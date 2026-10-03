#include <rainbow/patch_optics.hpp>
#include <rainbow/patch_accel.hpp>
#include <rainbow/patch_query.hpp>
#include <rainbow/cuda_error.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <stdexcept>

namespace rainbow
{
PatchOptics::PatchOptics(const CudaContext& cuda_context) noexcept
    : cuda_context_(cuda_context), module_(cuda_context), results_(cuda_context)
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

    cuda_context_.make_current();
    // Complete any use on this stream before reallocating our output. Other
    // streams are the caller's responsibility, as in the existing GPU classes.
    synchronize();
    has_result_ = false;
    results_.close();
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
        void* arguments[] = {&params};
        constexpr unsigned int threads = 128;
        const unsigned int blocks = n / threads + (n % threads != 0u ? 1u : 0u);
        has_pending_work_ = true;
        RAINBOW_CUDA_CHECK(cuLaunchKernel(function_, blocks, 1, 1, threads, 1, 1,
                                         0, cuda_context_.stream(), arguments, nullptr));
        synchronize();
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
