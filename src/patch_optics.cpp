#include <rainbow/patch_optics.hpp>
#include <rainbow/patch_accel.hpp>
#include <rainbow/patch_query.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/focal_phase.hpp>
#include <rainbow/rainbow_diffraction.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <stdexcept>

namespace rainbow
{
PatchOptics::PatchOptics(const CudaContext& cuda_context) noexcept
    : cuda_context_(cuda_context), module_(cuda_context), results_(cuda_context),
      focal_results_(cuda_context), transitions_(cuda_context), diffraction_results_(cuda_context)
{
}
PatchOptics::~PatchOptics() noexcept { static_cast<void>(close_noexcept()); }

void PatchOptics::load_module(const std::filesystem::path& fatbin_path)
{
    if(is_closed_) throw std::logic_error("PatchOptics is closed.");
    module_.load_fatbin(fatbin_path);
    function_ = module_.find_function("evaluate_patch_optics");
}
void PatchOptics::synchronize()
{
    cuda_context_.make_current();
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(cuda_context_.stream()));
    has_pending_work_ = false;
}

void PatchOptics::evaluate(const PatchAccel& source, const PatchQuery& query,
                          const RaindropTraceConfig& config)
{ evaluate_impl(source, query, config, nullptr); }

void PatchOptics::evaluate_wave(const PatchAccel& source, const PatchQuery& query,
                               const RaindropTraceConfig& config, const WaveOpticsSettings& settings)
{ evaluate_impl(source, query, config, &settings); }

void PatchOptics::evaluate_impl(const PatchAccel& source, const PatchQuery& query,
                               const RaindropTraceConfig& config, const WaveOpticsSettings* wave)
{
    if(is_closed_ || !function_) throw std::logic_error("Load the optical CUDA module first.");
    if(!source.has_result() || !query.has_result() || !query.matches_source(source))
        throw std::invalid_argument("Optics needs matching completed patch/query results.");
    PatchBuildLayout layout{};
    if(!PatchBuildLayout::try_make(config, layout)
       || layout.grid_width != source.layout().grid_width
       || layout.grid_height != source.layout().grid_height
       || layout.incident_area_drop2 != source.layout().incident_area_drop2)
        throw std::invalid_argument("Optics config does not match the source layout.");
    if(!(config.radius_mm > 0.0f && std::isfinite(config.radius_mm)
         && config.wavelength_nm > 0.0f && std::isfinite(config.wavelength_nm)))
        throw std::invalid_argument("Invalid optical unit metadata.");
    const auto n = query.statistics().directions;
    if(query.directions().element_count() != n || query.summaries().element_count() != n
       || query.offsets().element_count() != std::size_t(n) + 1u
       || query.host_directions().size() != n)
        throw std::invalid_argument("Query array lengths do not match.");

    DiffractionConfig new_diffraction{};
    FocalPhaseConfig focal_config{config.grid_width, config.grid_height, config.grid_half_extent, 0, {0,0,0,0}};
    if(wave)
    {
        if(!query.direction_grid() || n == 0u) throw std::invalid_argument("Wave optics requires a nonempty angular grid.");
        for(unsigned i=0;i<4;++i) focal_config.quarter_turn_offsets[i]=wave->focal_quarter_turn_offsets[i];
        if(!FocalPhase::layout_valid(focal_config)) throw std::invalid_argument("Invalid focal-phase layout/offsets.");
        double sigma=wave->primary_sigma_degrees;
        if(!(sigma>=0.0 && std::isfinite(sigma))) throw std::invalid_argument("Invalid diffraction sigma.");
        if(sigma==0.0 && !RainbowDiffraction::table_sigma_degrees(config.radius_mm, sigma))
            throw std::invalid_argument("Table II covers radii 0.1..1.0 mm only; specify --diffraction-sigma-deg for an explicit extension.");
        new_diffraction.primary_sigma_rad=sigma*RainbowDiffraction::pi/180.0;
        new_diffraction.transition_contrast=wave->transition_contrast;
        DiffractionParams check{};
        check.theta_count=query.direction_grid()->theta_count;check.phi_count=query.direction_grid()->phi_count;
        check.config=new_diffraction;
        if(!RainbowDiffraction::valid_config(check)) throw std::invalid_argument("Invalid diffraction/grid configuration.");
    }
    cuda_context_.make_current();
    // Complete any use on this stream before reallocating our output. Other
    // streams are the caller's responsibility, as in the existing GPU classes.
    synchronize();
    has_result_ = false; has_wave_result_ = false;
    results_.close(); focal_results_.close(); transitions_.close(); diffraction_results_.close();
    host_focal_.clear(); host_transitions_.clear(); host_diffraction_.clear(); wave_statistics_={};
    if(wave)
    {
        wave_settings_=*wave;diffraction_config_=new_diffraction;
        focal_results_.allocate(n);transitions_.allocate(n);diffraction_results_.allocate(n);
        host_focal_.resize(n);host_transitions_.resize(n);host_diffraction_.resize(n);
    }
    host_results_.assign(n, PatchOpticalResult{});
    host_directions_ = query.host_directions();
    config_ = config;
    has_grid_ = query.direction_grid() != nullptr;
    if(has_grid_) grid_ = *query.direction_grid();
    statistics_ = {};
    statistics_.directions = n;
    statistics_.missing_source_cells = query.statistics().missing_source_cells;
    statistics_.no_outgoing_source_cells = query.statistics().no_outgoing_source_cells;
    results_.allocate(n);
    if(n != 0u)
    {
        PatchOpticsParams params{};
        params.vertices = reinterpret_cast<const OutgoingVertex*>(source.source_vertices_address());
        params.patches = reinterpret_cast<const OutgoingPatch*>(source.patches().address());
        params.directions = reinterpret_cast<const Vec3*>(query.directions().address());
        params.offsets = reinterpret_cast<const std::uint64_t*>(query.offsets().address());
        params.summaries = reinterpret_cast<const PatchQuerySummary*>(query.summaries().address());
        params.hits = reinterpret_cast<const PatchQueryHit*>(query.hits().address());
        params.results = results_.data();
        params.hit_storage_count = query.hits().element_count();
        params.vertex_count = source.layout().vertices_per_path * 4u;
        params.patch_count = source.statistics().patch_count;
        params.direction_count = n;
        params.incident_direction = config.incident_direction;
        params.incident_basis_x = config.incident_basis_x;
        // The kernel takes one struct BY VALUE, not a pointer to host memory.
        constexpr unsigned int threads = 128;
        const unsigned int blocks = n / threads + (n % threads != 0u ? 1u : 0u);
        if(wave)
        {
            WaveOpticsParams wp{params, focal_config, focal_results_.data()};
            void* arguments[] = {&wp};
            const auto focal_function=module_.find_function("evaluate_focal_optics");
            has_pending_work_=true;
            RAINBOW_CUDA_CHECK(cuLaunchKernel(focal_function,blocks,1,1,threads,1,1,0,
                                             cuda_context_.stream(),arguments,nullptr));
            run_diffraction(); // same stream; no intermediate readback or new OptiX pipeline
        }
        else
        {
            void* arguments[] = {&params};
            has_pending_work_ = true;
            RAINBOW_CUDA_CHECK(cuLaunchKernel(function_, blocks, 1, 1, threads, 1, 1,
                                             0, cuda_context_.stream(), arguments, nullptr));
        }
        synchronize();
        if(wave)
        {
            focal_results_.download(host_focal_);transitions_.download(host_transitions_);
            diffraction_results_.download(host_diffraction_);
        }
        results_.download(std::span<PatchOpticalResult>{host_results_});
    }
    for(const auto& r : host_results_)
    {
        statistics_.known_hits_complete_directions += r.known_hits_complete();
        statistics_.pending_directions += (r.flags & patch_optical_pending_mask) != 0u;
        statistics_.error_directions += (r.flags & patch_optical_error_mask) != 0u;
        statistics_.evaluated_hits += r.evaluated_hits;
        statistics_.rejected_hits += r.rejected_hits;
    }
    // A true result means that evaluation finished, NOT that all physics or all
    // source cells are complete. CLI inspects errors only after writing the CSV.
    has_result_ = true;
    if(wave)
    {
        for(std::size_t i=0;i<host_focal_.size();++i)
        {
            const auto& f=host_focal_[i];const auto& d=host_diffraction_[i];
            const auto t=host_transitions_[i].flags;
            wave_statistics_.focal_errors+=(f.flags&focal_error_mask)!=0u;
            wave_statistics_.focal_pending+=(f.flags&focal_pending_mask)!=0u;
            wave_statistics_.corrected_hits+=f.corrected_hits;
            wave_statistics_.extra_quarter_turn_hits+=f.extra_quarter_turn_hits;
            wave_statistics_.primary_transitions+=(t&TransitionPrimary)!=0u;
            wave_statistics_.secondary_transitions+=(t&TransitionSecondary)!=0u;
            wave_statistics_.unknown_transitions+=(t&TransitionUnknown)!=0u;
            wave_statistics_.filtered_directions+=(d.flags&DiffractionFiltered)!=0u;
            wave_statistics_.unavailable_directions+=!d.valid();
            wave_statistics_.underresolved_directions+=(d.flags&DiffractionUnderresolved)!=0u;
        }
        has_wave_result_=true;
    }
}

void PatchOptics::run_diffraction()
{
    DiffractionParams p{};
    p.focal=focal_results_.data();p.transitions=transitions_.data();p.results=diffraction_results_.data();
    p.theta_count=grid_.theta_count;p.phi_count=grid_.phi_count;p.config=diffraction_config_;
    void* args[]={&p};
    constexpr unsigned threads=128;
    const auto n=p.theta_count*p.phi_count;
    const unsigned blocks=n/threads+(n%threads!=0u?1u:0u);
    const auto detect=module_.find_function("detect_rainbow_transitions");
    const auto filter=module_.find_function("filter_rainbow_diffraction");
    has_pending_work_=true;
    RAINBOW_CUDA_CHECK(cuLaunchKernel(detect,blocks,1,1,threads,1,1,0,cuda_context_.stream(),args,nullptr));
    RAINBOW_CUDA_CHECK(cuLaunchKernel(filter,blocks,1,1,threads,1,1,0,cuda_context_.stream(),args,nullptr));
}

void PatchOptics::write_csv(const std::filesystem::path& path) const
{
    if(!has_result_) throw std::logic_error("No optical result to save.");
    std::ofstream out(path);
    if(!out) throw std::runtime_error("Cannot open optical CSV.");
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "# format=rainbow_patch_optics_v1\n"
        << "# quantity=regular_partial_model_angular_density_NOT_phase_function\n"
        << "# optical_complete=false\n# source_coverage_certified=false\n"
        << "# focal_line_phase_applied=false\n# diffraction_applied=false\n# normalized=false\n"
        << "# input_states=one_coherent_Jones_state\n# phasor_convention=exp(+i*2*pi*q)\n"
        << "# interpolation=transport_then_bilinear_field_and_path_then_propagation_phase\n"
        << "# transport=shortest_great_circle_minimum_rotation_explicit_project_convention\n"
        << "# field_components=s_perpendicular,p_outgoing_cross_s\n"
        << "# density_units=input_field_squared*drop_unit_squared_per_sr\n"
        << "# physical_area_factor_mm2=" << double(config_.radius_mm) * config_.radius_mm << '\n'
        << "# rejected_hit_policy=retain_regular_partial_sum_and_mark_incomplete\n"
        << "# numerical_error_policy=all_optical_values_nan\n"
        << "# partial_coherent_intensity_is_NOT_a_lower_bound=true\n"
        << "# missing_source_cells=" << statistics_.missing_source_cells
        << "\n# no_outgoing_source_cells=" << statistics_.no_outgoing_source_cells << '\n'
        << "# radius_mm=" << config_.radius_mm << "\n# wavelength_nm=" << config_.wavelength_nm
        << "\n# exterior_index=" << config_.exterior_index << "\n# interior_index=" << config_.interior_index
        << "\n# incident_grid=" << config_.grid_width << ',' << config_.grid_height
        << "\n# incident_grid_half_extent_drop=" << config_.grid_half_extent
        << "\n# reference_distance_drop=" << config_.reference_distance
        << "\n# outgoing_reference_distance_drop=" << config_.outgoing_reference_distance << '\n'
        << "# incident_direction=" << config_.incident_direction.x << ',' << config_.incident_direction.y << ',' << config_.incident_direction.z
        << "\n# incident_basis_x=" << config_.incident_basis_x.x << ',' << config_.incident_basis_x.y << ',' << config_.incident_basis_x.z
        << "\n# incident_field=" << config_.incident_field.x.real << ',' << config_.incident_field.x.imag << ','
        << config_.incident_field.y.real << ',' << config_.incident_field.y.imag << '\n';
    out << "# coefficients=";
    for(unsigned i = 0; i < 8; ++i) out << (i ? "," : "") << config_.shape.coefficients[i];
    out << '\n';
    if(has_grid_)
        out << "# theta_count=" << grid_.theta_count << "\n# phi_count=" << grid_.phi_count
            << "\n# order=theta_major_phi_minor\n# angle_units=radians\n"
            << "# theta=acos(incident_propagation_direction_dot_outgoing_direction)\n"
            << "# phi=atan2(outgoing_dot_incident_basis_y,outgoing_dot_incident_basis_x)\n";
    out << "direction_id,theta_index,phi_index,theta_rad,phi_rad,solid_angle_sr,wx,wy,wz,"
           "known_hits_complete,hit_count,evaluated_hits,rejected_hits,refinement_hits,boundary_hits,singular_hits,"
           "flags,query_flags,first_problem_patch_id,"
           "regular_partial_incoherent_s,regular_partial_incoherent_p,regular_partial_incoherent_total,"
           "regular_partial_path_s_real,regular_partial_path_s_imag,regular_partial_path_p_real,regular_partial_path_p_imag,"
           "regular_partial_path_s,regular_partial_path_p,regular_partial_path_total\n";
    for(std::size_t i = 0; i < host_results_.size(); ++i)
    {
        const auto& r = host_results_[i];
        const auto& f = r.regular_partial_path_field;
        const auto& w = host_directions_[i];
        out << i << ',';
        if(has_grid_)
        {
            const auto row = static_cast<std::uint32_t>(i / grid_.phi_count);
            const auto col = static_cast<std::uint32_t>(i % grid_.phi_count);
            out << row << ',' << col << ',' << grid_.theta(row) << ',' << grid_.phi(col) << ',' << grid_.solid_angle(row) << ',';
        }
        else out << "-1,-1,nan,nan,nan,";
        out << w.x << ',' << w.y << ',' << w.z << ',' << r.known_hits_complete() << ','
            << r.hit_count << ',' << r.evaluated_hits << ',' << r.rejected_hits << ','
            << r.refinement_hits << ',' << r.boundary_hits << ',' << r.singular_hits << ','
            << r.flags << ',' << r.query_flags << ',';
        if(r.first_problem_patch_id == 0xffffffffu) out << -1; else out << r.first_problem_patch_id;
        out << ',' << r.regular_partial_incoherent_s << ',' << r.regular_partial_incoherent_p << ','
            << r.regular_partial_incoherent_s + r.regular_partial_incoherent_p << ','
            << f.s_real << ',' << f.s_imag << ',' << f.p_real << ',' << f.p_imag << ','
            << r.regular_partial_path_s << ',' << r.regular_partial_path_p << ','
            << r.regular_partial_path_s + r.regular_partial_path_p << '\n';
    }
    out.flush();
    if(!out) throw std::runtime_error("Writing optical CSV failed.");
}

void PatchOptics::write_wave_csv(const std::filesystem::path& path) const
{
    if(!has_wave_result_ || !has_result_ || !has_grid_) throw std::logic_error("No wave-optics grid result.");
    std::ofstream out(path);
    if(!out) throw std::runtime_error("Cannot open wave-optics CSV.");
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "# format=rainbow_wave_optics_v1\n"
        << "# quantity=model_angular_density_NOT_certified_phase_function\n"
        << "# optical_effects_connected=true\n# optical_complete=false\n# source_coverage_certified=false\n"
        << "# normalized=false\n# focal_line_phase_applied=true\n# diffraction_stage_executed=true\n"
        << "# diffraction_applied_locally=true\n# diffraction_width_wavelength_scaling=none\n"
        << "# phasor_convention=exp(+i*2*pi*q)\n# focal_rule=radial_cell_center_sign_plus_family_offsets\n"
        << "# focal_absolute_family_offsets=explicit_project_convention_NOT_specified_in_paper\n"
        << "# focal_quarter_turn_offsets=";
    for(unsigned f=0;f<4;++f)out<<(f?",":"")<<wave_settings_.focal_quarter_turn_offsets[f];
    out << "\n# focal_derivative=radial_derivative_of_bilinear_theta_at_emitting_cell_center\n"
        << "# focal_invalid_policy=nan_no_complete_looking_partial_sum\n"
        << "# diffraction_model=project_local_spherical_gaussian_v1\n"
        << "# diffraction_not_author_code_identical=true\n"
        << "# diffraction_sigma_source=" << (wave_settings_.primary_sigma_degrees>0?"explicit_user_override":"Table_II_linear_radius_interpolation")
        << "\n# diffraction_primary_sigma_deg="<<diffraction_config_.primary_sigma_rad*180/RainbowDiffraction::pi
        << "\n# diffraction_secondary_sigma_deg="<<2*diffraction_config_.primary_sigma_rad*180/RainbowDiffraction::pi
        << "\n# diffraction_support_sigma="<<diffraction_config_.support_sigma
        << "\n# diffraction_inner_blend_sigma="<<diffraction_config_.inner_blend_sigma
        << "\n# diffraction_transition_contrast="<<diffraction_config_.transition_contrast
        << "\n# diffraction_detection=sharp_total_intensity_jump_AND_primary_enter_or_secondary_leave_support\n"
        << "# diffraction_kernel=exp(-angular_distance_squared/(2*sigma_squared))*cell_solid_angle\n"
        << "# diffraction_components=local_s_p_scalar_intensities_NOT_complex_field_filter\n"
        << "# diffraction_holes=propagate_unavailable_do_NOT_renormalize_over_valid_only\n"
        << "# global_energy_conservation_guaranteed=false\n# independent_peak_or_energy_fitting=false\n"
        << "# angular_quadrature=point_samples_times_exact_cell_solid_angles\n"
        << "# density_units=input_field_squared*drop_unit_squared_per_sr\n"
        << "# physical_area_factor_mm2="<<double(config_.radius_mm)*config_.radius_mm
        << "\n# radius_mm="<<config_.radius_mm<<"\n# wavelength_nm="<<config_.wavelength_nm
        << "\n# exterior_index="<<config_.exterior_index<<"\n# interior_index="<<config_.interior_index
        << "\n# incident_grid="<<config_.grid_width<<','<<config_.grid_height
        << "\n# incident_grid_half_extent_drop="<<config_.grid_half_extent
        << "\n# reference_distance_drop="<<config_.reference_distance
        << "\n# outgoing_reference_distance_drop="<<config_.outgoing_reference_distance
        << "\n# incident_direction="<<config_.incident_direction.x<<','<<config_.incident_direction.y<<','<<config_.incident_direction.z
        << "\n# incident_basis_x="<<config_.incident_basis_x.x<<','<<config_.incident_basis_x.y<<','<<config_.incident_basis_x.z
        << "\n# incident_field="<<config_.incident_field.x.real<<','<<config_.incident_field.x.imag<<','
        <<config_.incident_field.y.real<<','<<config_.incident_field.y.imag
        << "\n# coefficients=";
    for(unsigned f=0;f<8;++f)out<<(f?",":"")<<config_.shape.coefficients[f];
    out << "\n# missing_source_cells="<<statistics_.missing_source_cells
        << "\n# no_outgoing_source_cells="<<statistics_.no_outgoing_source_cells
        << "\n# theta_count="<<grid_.theta_count<<"\n# phi_count="<<grid_.phi_count
        << "\n# order=theta_major_phi_minor\n# angle_units=radians\n"
        << "# theta=acos(incident_propagation_direction_dot_outgoing_direction)\n"
        << "# phi=atan2(outgoing_dot_incident_basis_y,outgoing_dot_incident_basis_x)\n"
        << "# focal_flag_bits=1:input_error,2:input_pending,4:derivative_pending,8:geometry_error,16:arithmetic_error\n"
        << "# diffraction_flag_bits=1:filtered,2:underresolved,4:input_unavailable,8:stencil_unavailable,16:detection_unavailable,32:invalid_data\n"
        << "direction_id,theta_index,phi_index,theta_rad,phi_rad,solid_angle_sr,wx,wy,wz,"
        << "known_hits_complete,hit_count,evaluated_hits,rejected_hits,optical_flags,"
        << "focal_valid,focal_flags,focal_corrected_hits,focal_extra_quarter_turn_hits,focal_first_problem_patch_id,"
        << "R_hits,TT_hits,TRT_hits,TRRT_hits,transition_flags,diffraction_valid,diffraction_flags,"
        << "diffraction_transition_kind,diffraction_sigma_rad,diffraction_blend,"
        << "incoherent_s,incoherent_p,incoherent_total,path_s,path_p,path_total,"
        << "focal_s_real,focal_s_imag,focal_p_real,focal_p_imag,focal_s,focal_p,focal_total,"
        << "diffraction_s,diffraction_p,diffraction_total\n";
    for(std::size_t i=0;i<host_results_.size();++i)
    {
        const auto& r=host_results_[i];const auto& f=host_focal_[i];const auto& d=host_diffraction_[i];
        const auto& w=host_directions_[i];
        const auto row=static_cast<std::uint32_t>(i/grid_.phi_count),col=static_cast<std::uint32_t>(i%grid_.phi_count);
        out<<i<<','<<row<<','<<col<<','<<grid_.theta(row)<<','<<grid_.phi(col)<<','<<grid_.solid_angle(row)<<','
           <<w.x<<','<<w.y<<','<<w.z<<','<<r.known_hits_complete()<<','<<r.hit_count<<','<<r.evaluated_hits<<','
           <<r.rejected_hits<<','<<r.flags<<','<<f.valid()<<','<<f.flags<<','<<f.corrected_hits<<','
           <<f.extra_quarter_turn_hits<<',';
        if(f.first_problem_patch_id==0xffffffffu)out<<-1;else out<<f.first_problem_patch_id;
        for(unsigned j=0;j<4;++j)out<<','<<f.family_hits[j];
        out<<','<<host_transitions_[i].flags<<','<<d.valid()<<','<<d.flags<<','<<d.transition_kind<<','<<d.sigma_rad<<','<<d.blend
           <<','<<r.regular_partial_incoherent_s<<','<<r.regular_partial_incoherent_p<<','
           <<r.regular_partial_incoherent_s+r.regular_partial_incoherent_p
           <<','<<r.regular_partial_path_s<<','<<r.regular_partial_path_p<<','<<r.regular_partial_path_s+r.regular_partial_path_p
           <<','<<f.field.s_real<<','<<f.field.s_imag<<','<<f.field.p_real<<','<<f.field.p_imag
           <<','<<f.intensity_s<<','<<f.intensity_p<<','<<f.intensity_s+f.intensity_p
           <<','<<d.intensity_s<<','<<d.intensity_p<<','<<d.intensity_s+d.intensity_p<<'\n';
    }
    out.flush();if(!out)throw std::runtime_error("Writing wave-optics CSV failed.");
}

bool PatchOptics::close_noexcept() noexcept
{
    if(is_closed_) return true;
    if(!detail::report_cuda_cleanup_result(cuCtxSetCurrent(cuda_context_.handle()),
                                           "cuCtxSetCurrent(optics cleanup)")) return false;
    bool ok = true;
    if(has_pending_work_)
    {
        ok = detail::report_cuda_cleanup_result(cuStreamSynchronize(cuda_context_.stream()),
                                               "cuStreamSynchronize(optics cleanup)") && ok;
        has_pending_work_ = false;
    }
    ok = diffraction_results_.close_noexcept() && ok;
    ok = transitions_.close_noexcept() && ok;
    ok = focal_results_.close_noexcept() && ok;
    has_wave_result_ = false;
    ok = results_.close_noexcept() && ok;
    ok = module_.close_noexcept() && ok;
    function_ = nullptr; has_result_ = false; is_closed_ = true;
    return ok;
}
void PatchOptics::close()
{
    if(!close_noexcept()) throw std::runtime_error("Patch optics cleanup failed; see stderr.");
}
} // namespace rainbow
