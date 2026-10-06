#include <rainbow/patch_failure_report.hpp>
#include <bit>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace rainbow
{
namespace
{
static_assert(sizeof(float)==4 && std::numeric_limits<float>::is_iec559);
static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits==53);
void number(std::ostream& out,const double value)
{if(std::isfinite(value))out<<value;else out<<"null";}
void bits(std::ostream& out,const float value) {out<<std::bit_cast<std::uint32_t>(value);}
void vec_bits(std::ostream& out,const Vec3 p)
{out<<'[';bits(out,p.x);out<<',';bits(out,p.y);out<<',';bits(out,p.z);out<<']';}
void field_bits(std::ostream& out,const Field32 f)
{out<<'[';bits(out,f.x.real);out<<',';bits(out,f.x.imag);out<<',';bits(out,f.y.real);out<<',';bits(out,f.y.imag);out<<']';}
void field(std::ostream& out,const OpticalField64 f)
{out<<'[';number(out,f.s_real);out<<',';number(out,f.s_imag);out<<',';number(out,f.p_real);out<<',';number(out,f.p_imag);out<<']';}
void bounds(std::ostream& out,const std::array<PatchRegularityAudit::Bound,4>& x)
{
    out<<'[';
    for(unsigned k=0;k<4;++k){if(k)out<<',';out<<'[';number(out,x[k].lo);out<<',';number(out,x[k].hi);out<<']';}
    out<<']';
}
void vertex(std::ostream& out,const OutgoingVertex& v,const std::uint32_t id)
{
    out<<"{\"vertex_index\":"<<id<<",\"position_bits\":";vec_bits(out,v.position_drop);
    out<<",\"direction_bits\":";vec_bits(out,v.direction_drop);
    out<<",\"basis_x_bits\":";vec_bits(out,v.basis_x);
    out<<",\"field_bits\":";field_bits(out,v.field);
    out<<",\"turns\":"<<v.optical_cycles.turns<<",\"fraction_bits\":";bits(out,v.optical_cycles.fraction);
    out<<",\"status\":"<<static_cast<unsigned>(v.status)<<",\"diagnostics\":"<<v.diagnostics<<'}';
}
void inspect_patch(std::ostream& out,const FailedPatchWitness& w,const RaindropTraceConfig& config)
{
    const auto& p=w.patch;
    const auto nx=config.grid_width-1u,ny=config.grid_height-1u;
    const auto cells=std::uint64_t(nx)*ny;
    const auto cell=p.patch_id%cells;
    const double U=(2*(double(cell%nx)+0.5)/nx-1)*config.grid_half_extent;
    const double V=(2*(double(cell/nx)+0.5)/ny-1)*config.grid_half_extent;
    out<<"{\"compact_index\":"<<w.compact_index<<",\"patch_id\":"<<p.patch_id
       <<",\"family\":"<<p.patch_id/cells<<",\"cell_x\":"<<cell%nx<<",\"cell_y\":"<<cell/nx
       <<",\"stored_status\":"<<static_cast<unsigned>(p.status)<<",\"area_bits\":";
    bits(out,p.incident_area_drop2);out<<",\"stored_omega_bits\":";bits(out,p.signed_solid_angle_sr);
    out<<",\"cell_center\":[";number(out,U);out<<',';number(out,V);out<<"]";
    out<<",\"vertices\":[";
    BilinearPatchGeometry g{};
    for(unsigned k=0;k<4;++k)
    {if(k)out<<',';vertex(out,w.vertices[k],p.vertex_indices[k]);g.corners[k]=w.vertices[k].direction_drop;}
    const auto a=PatchRegularityAudit::inspect(g);
    out<<"],\"host_audit\":{\"finite_input\":"<<(a.finite_input?"true":"false")
       <<",\"orientation32\":"<<a.host_orientation32<<",\"interval_orientation64\":"<<a.interval_orientation64
       <<",\"hemisphere_axis_bits\":";vec_bits(out,a.hemisphere_axis);
    out<<",\"jacobian32\":";bounds(out,a.jacobian32);
    out<<",\"jacobian64\":";bounds(out,a.jacobian64);
    out<<",\"hemisphere32\":";bounds(out,a.hemisphere32);
    out<<",\"hemisphere64\":";bounds(out,a.hemisphere64);
    out<<",\"jacobian_value64\":[";
    for(unsigned k=0;k<4;++k){if(k)out<<',';number(out,a.jacobian_value64[k]);}
    out<<"],\"signed_omega32_bits\":";bits(out,a.host_signed_omega32);
    out<<",\"signed_omega64\":";number(out,a.signed_omega64);
    out<<"}}";
}
}

PatchFailureReport::PatchFailureReport(PatchFailureSnapshot snapshot):snapshot_(std::move(snapshot)) {}

void PatchFailureReport::write_json(const std::filesystem::path& path) const
{
    const auto& s=snapshot_;
    const auto& c=s.config;
    if(c.grid_width<2 || c.grid_height<2 || s.theta_count==0 || s.phi_count==0)
        throw std::invalid_argument("Patch report requires nonempty valid grids.");
    std::map<std::uint32_t,const FailedPatchWitness*> patch_by_compact;
    std::set<std::uint32_t> patch_ids,direction_ids;
    for(const auto& p:s.patches)
    {
        if(!patch_by_compact.emplace(p.compact_index,&p).second || !patch_ids.insert(p.patch.patch_id).second)
            throw std::invalid_argument("Duplicate patch witness.");
    }
    const auto count=std::uint64_t(s.theta_count)*s.phi_count;
    for(const auto& d:s.directions)
    {
        if(d.direction_id>=count || !direction_ids.insert(d.direction_id).second || d.hits.size()!=d.query.hit_count)
            throw std::invalid_argument("Invalid/duplicate direction witness or hit count.");
        for(const auto& h:d.hits)
        {
            const auto it=patch_by_compact.find(h.compact_index);
            if(it==patch_by_compact.end() || it->second->patch.patch_id!=h.patch_id)
                throw std::invalid_argument("Hit refers to an absent/different witness.");
        }
    }
    if(path.empty() || std::filesystem::exists(path))
        throw std::invalid_argument("Patch report output must be a NEW file; preserve the previous report or choose another name.");
    auto temporary=path;temporary+=".part";
    if(std::filesystem::exists(temporary))
        throw std::invalid_argument("A previous .part report exists; inspect/remove it or choose another report path.");
    std::ofstream out(temporary,std::ios::binary);
    if(!out)throw std::runtime_error("Cannot create patch failure report.");
    out.imbue(std::locale::classic());out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    const char* origin=s.origin==PatchWitnessOrigin::GpuCapture?"gpu_buffer_snapshot":
                       (s.origin==PatchWitnessOrigin::CpuRetrace?"cpu_retrace_NOT_gpu":"synthetic_test");
    out<<"{\n\"format\":\"rainbow_patch_failure_witness_v1\",\n\"origin\":\""<<origin<<"\",\n"
       <<"\"integration_baseline_commit\":\"8c8deabef6fb10717308491b1e992a266b46e3d0\",\n"
       <<"\"changes_numerical_results\":false,\n"
       <<"\"float_encoding\":\"IEEE754_binary32_bits_as_uint32\",\n"
       <<"\"scope\":\"all_recorded_hits_of_failed_directions_plus_named_first_problem_patches\",\n"
       <<"\"limits\":\"No source coverage certificate; query failures without hits expose only recorded first_problem IDs; host replay is not a GPU predicate trace.\",\n"
       <<"\"config\":{\"grid_width\":"<<c.grid_width<<",\"grid_height\":"<<c.grid_height
       <<",\"theta_count\":"<<s.theta_count<<",\"phi_count\":"<<s.phi_count
       <<",\"radius_bits\":";bits(out,c.radius_mm);
    out<<",\"wavelength_bits\":";bits(out,c.wavelength_nm);
    out<<",\"ior_bits\":[";bits(out,c.exterior_index);out<<',';bits(out,c.interior_index);
    out<<"],\"half_extent_bits\":";bits(out,c.grid_half_extent);
    out<<",\"reference_distance_bits\":[";bits(out,c.reference_distance);out<<',';bits(out,c.outgoing_reference_distance);
    out<<"],\"incident_direction_bits\":";vec_bits(out,c.incident_direction);
    out<<",\"incident_basis_x_bits\":";vec_bits(out,c.incident_basis_x);
    out<<",\"incident_field_bits\":";field_bits(out,c.incident_field);
    out<<",\"shape_coefficients_bits\":[";
    for(unsigned k=0;k<8;++k){if(k)out<<',';bits(out,c.shape.coefficients[k]);}
    out<<"],\"focal_offsets\":[";
    for(unsigned k=0;k<4;++k){if(k)out<<',';out<<s.wave_settings.focal_quarter_turn_offsets[k];}
    out<<"]},\n\"summary\":{\"total_directions\":"<<count
       <<",\"selected_directions\":"<<s.directions.size()<<",\"selected_patches\":"<<s.patches.size()
       <<",\"stored_patch_count\":"<<s.stored_patch_count
       <<",\"missing_source_cells\":"<<s.missing_source_cells
       <<",\"no_outgoing_source_cells\":"<<s.no_outgoing_source_cells
       <<",\"query_error_directions\":"<<s.query_error_directions
       <<",\"optical_incomplete_directions\":"<<s.optical_incomplete_directions
       <<",\"focal_unavailable_directions\":"<<s.focal_unavailable_directions
       <<",\"diffraction_unavailable_directions\":"<<s.diffraction_unavailable_directions
       <<",\"nonfinite_directions\":"<<s.nonfinite_directions
       <<",\"underresolved_directions\":"<<s.underresolved_directions
       <<",\"gpu_read_bytes\":"<<s.gpu_read_bytes<<",\"gpu_read_calls\":"<<s.gpu_read_calls<<"},\n";
    out<<"\"patches\":[\n";
    for(std::size_t i=0;i<s.patches.size();++i){if(i)out<<",\n";inspect_patch(out,s.patches[i],c);}
    out<<"\n],\n\"directions\":[\n";
    constexpr double pi=3.141592653589793238462643383279502884;
    for(std::size_t i=0;i<s.directions.size();++i)
    {
        if(i)out<<",\n";
        const auto& d=s.directions[i];const auto& q=d.query;const auto& o=d.optical;
        const auto& f=d.focal;const auto& df=d.diffraction;
        const auto row=d.direction_id/s.phi_count,col=d.direction_id%s.phi_count;
        out<<"{\"direction_id\":"<<d.direction_id<<",\"theta_rad\":"<<pi*(double(row)+0.5)/s.theta_count
           <<",\"phi_rad\":"<<-pi+2*pi*(double(col)+0.5)/s.phi_count<<",\"direction_bits\":";vec_bits(out,d.direction);
        out<<",\"query_flags\":"<<q.flags<<",\"query_first_problem\":"<<q.first_problem_patch_id
           <<",\"candidate_count\":"<<q.candidate_count<<",\"hit_count\":"<<q.hit_count
           <<",\"optical_flags\":"<<o.flags<<",\"optical_first_problem\":"<<o.first_problem_patch_id
           <<",\"evaluated_hits\":"<<o.evaluated_hits<<",\"rejected_hits\":"<<o.rejected_hits
           <<",\"focal_flags\":"<<f.flags<<",\"focal_first_problem\":"<<f.first_problem_patch_id
           <<",\"focal_corrected_hits\":"<<f.corrected_hits
           <<",\"focal_extra_quarter_turn_hits\":"<<f.extra_quarter_turn_hits
           <<",\"family_hits\":[";
        for(unsigned k=0;k<4;++k){if(k)out<<',';out<<f.family_hits[k];}
        out<<"],\"diffraction_flags\":"<<df.flags<<",\"diffraction_transition_kind\":"<<df.transition_kind
           <<",\"diffraction_sigma_rad\":";number(out,df.sigma_rad);
        out<<",\"diffraction_blend\":";number(out,df.blend);
        out<<",\"incoherent\":[";number(out,o.regular_partial_incoherent_s);out<<',';number(out,o.regular_partial_incoherent_p);
        out<<"],\"path_field\":";field(out,o.regular_partial_path_field);
        out<<",\"path_intensity\":[";number(out,o.regular_partial_path_s);out<<',';number(out,o.regular_partial_path_p);
        out<<"],\"focal_field\":";field(out,f.field);
        out<<",\"focal_intensity\":[";number(out,f.intensity_s);out<<',';number(out,f.intensity_p);
        out<<"],\"diffraction_intensity\":[";number(out,df.intensity_s);out<<',';number(out,df.intensity_p);
        out<<"],\"hits\":[";
        for(std::size_t k=0;k<d.hits.size();++k)
        {
            if(k)out<<',';
            const auto& h=d.hits[k];
            const auto* w=patch_by_compact.at(h.compact_index);
            BilinearPatchGeometry g{};for(unsigned j=0;j<4;++j)g.corners[j]=w->vertices[j].direction_drop;
            out<<"{\"compact_index\":"<<h.compact_index<<",\"patch_id\":"<<h.patch_id<<",\"root_index\":"<<h.root_index
               <<",\"flags\":"<<h.flags<<",\"uvtr_bits\":[";
            bits(out,h.u);out<<',';bits(out,h.v);out<<',';bits(out,h.t);out<<',';bits(out,h.residual);
            out<<"],\"jacobian_at_hit64\":";number(out,PatchRegularityAudit::jacobian_at(g,h.u,h.v));out<<'}';
        }
        out<<"]}";
    }
    out<<"\n],\n\"absent_problem_patch_ids\":[";
    for(std::size_t k=0;k<s.absent_problem_patch_ids.size();++k){if(k)out<<',';out<<s.absent_problem_patch_ids[k];}
    out<<"],\n\"report_complete\":true\n}\n";
    out.flush();if(!out)throw std::runtime_error("Writing patch failure report failed; incomplete data remains in .part.");
    out.close();if(!out)throw std::runtime_error("Closing patch failure report failed.");
    std::filesystem::rename(temporary,path);
}
} // namespace rainbow
