#include <ice_crystal/phase_generate.hpp>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <locale>
#include <stdexcept>

#ifndef ICE_PHASE_SOURCE_FINGERPRINT
#define ICE_PHASE_SOURCE_FINGERPRINT "unknown"
#endif
namespace iceCrystal
{
namespace
{
void quote(std::ostream& o,std::string_view s)
{
    o<<'"';constexpr char hex[]="0123456789abcdef";
    for(unsigned char c:s)
    {
        if(c=='"'||c=='\\')o<<'\\'<<char(c);
        else if(c<32)o<<"\\u00"<<hex[c>>4]<<hex[c&15];
        else o<<char(c);
    }
    o<<'"';
}
std::ofstream output(const std::filesystem::path& p)
{
    std::ofstream out(p,std::ios::binary|std::ios::trunc);
    if(!out)throw std::runtime_error("Cannot open metadata output.");
    out.exceptions(std::ios::badbit|std::ios::failbit);out.imbue(std::locale::classic());out<<std::setprecision(17);return out;
}
void vector(std::ostream& f,Vec3 v){f<<'['<<v.x<<','<<v.y<<','<<v.z<<']';}
}
void validate_generation_settings(const PhaseGenerationSettings& s)
{
    validate_phase_grid(s.grid);validate_cdf_policy(s.cdf_policy);
    static_cast<void>(make_trace_settings(s.prism,s.grid.ki,s.exterior_index,s.interior_index,1,0,s.max_internal_hits,s.tail_tolerance));
    if(!(s.wavelength_nm>0)||!finite_value(s.wavelength_nm)||!s.batch_size||!s.max_output_records||
       !s.out_nt||!s.out_np||s.out_nt>s.grid.nt||s.out_np>s.grid.np||s.grid.nt%s.out_nt||s.grid.np%s.out_np)
        throw std::invalid_argument("Invalid wavelength, capacity, or coarsening dimensions.");
    for(double tolerance:{s.max_unresolved_fraction,s.max_balance_error,s.max_histogram_relative_error})
        if(!finite_value(tolerance)||tolerance<0||tolerance>=1)throw std::invalid_argument("Invalid quality tolerance.");
}
TraceSettings node_trace_settings(const PhaseGenerationSettings& s,const OrientationNode& n)
{
    if(!n.rotation.valid()||!finite_value(n.weight)||!(n.weight>0)||!n.samples)throw std::invalid_argument("Invalid orientation node.");
    return make_trace_settings(s.prism,n.rotation.inverse(s.grid.ki),s.exterior_index,s.interior_index,
        n.samples,n.seed,s.max_internal_hits,s.tail_tolerance);
}
double node_area(const PhaseGenerationSettings& s,const OrientationNode& n)
{
    const auto c=node_trace_settings(s,n);return s.prism.projected_area_mm2(c.incident_direction);
}
TraceSummary summarize_phase_rays(std::span<const RayAudit> rays,double tolerance)
{
    auto result=summarize(rays,tolerance);
    for(const auto& r:rays)
        result.accepted=result.accepted&&finite_value(r.max_interface_balance_error)&&
            r.max_interface_balance_error>=0&&r.max_interface_balance_error<=tolerance;
    return result;
}
void record_batch_audit(CompositionAudit& a,const TraceSummary& t,double scale)
{
    merge_summary(a.trace,t);
    a.escaped_area.add(scale*t.escaped_power_sum);a.unresolved_area.add(scale*t.unresolved_power_sum);
    require_accepted(t); // a bad ray is never omitted while retaining a completed record
}
void validate_composition(CompositionAudit& a,const PhaseGenerationSettings& s,const OrientationPlan& plan)
{
    require_accepted(a.trace);
    if(a.trace.rays!=plan.total_samples||a.orientations_done!=plan.nodes.size()||
       a.histogram.samples!=a.trace.outputs||a.histogram.invalid_samples)throw std::runtime_error("Incomplete composition / output-count mismatch.");
    const double area=a.expected_area.value,escaped=a.escaped_area.value,unresolved=a.unresolved_area.value;
    if(!(area>0)||!finite_value(area)||!(escaped>0)||!finite_value(escaped)||!finite_value(unresolved)||unresolved<0)
        throw std::runtime_error("Invalid integrated cross sections.");
    if(unresolved/area>s.max_unresolved_fraction)throw std::runtime_error("Unresolved light exceeds explicit global budget.");
    if(std::fabs((escaped+unresolved)/area-1)>s.max_balance_error)throw std::runtime_error("Orientation-weighted energy balance failed.");
    if(!(a.grid_moments.mass>0)||!finite_value(a.grid_moments.mass)||!finite_value(a.grid_moments.axial)||
       !(a.histogram.point_mass>0)||!finite_value(a.histogram.point_mass)||!finite_value(a.histogram.point_axial))
        throw std::runtime_error("Nonfinite histogram moments.");
    a.histogram_relative_error=std::max(std::fabs(a.grid_moments.mass/escaped-1),std::fabs(a.histogram.point_mass/escaped-1));
    if(a.histogram_relative_error>s.max_histogram_relative_error)throw std::runtime_error("Histogram mass differs from independent escaped-ray audit.");
}
void write_phase_metadata(const std::filesystem::path& path,const PhaseGenerationSettings& s,
    const OrientationPlan& plan,const CompositionAudit& a,const CdfReport& r,const char* backend)
{
    auto f=output(path);const double area=a.expected_area.value;
    f<<"{\n  \"schema\": \"rainbow.phase_cdf.numpy.v2\",\n  \"complete\": true,\n"
     <<"  \"generator_version\": \"ice.hex_phase.m2.v1\",\n  \"backend\": ";quote(f,backend);
    f<<",\n  \"base_revision\": \"a9538941dcb9df2de6a4130fa66f625f0133c949 + hex-trace-m1\",\n"
     <<"  \"source_fingerprint_at_configure\": \""<<ICE_PHASE_SOURCE_FINGERPRINT<<"\",\n"
     <<"  \"dtype\": \"<f8\",\n  \"order\": \"C\",\n"
     <<"  \"theta_count\": "<<s.out_nt<<",\n  \"phi_count\": "<<s.out_np<<",\n"
     <<"  \"density_measure\": \"solid_angle_sr\",\n  \"cell_model\": \"constant_density_per_spherical_cell\",\n"
     <<"  \"coordinates\": \"u=(1-cos(theta))/2; v=(phi+pi)/(2*pi)\",\n"
     <<"  \"cdf_axis_order\": [\"phi_cell\",\"theta_edge\"],\n"
     <<"  \"normalization\": \"integral_p_domega_equals_one\",\n"
     <<"  \"direction_convention\": \"physical_propagation\",\n"
     <<"  \"coordinate_contract_id\": \"ice.phase_cdf.ensemble_coordinates.v1\",\n"
     <<"  \"sampling_frame_layout\": \"rows_xyz_columns_e0_e1_k\",\n"
     <<"  \"sampling_frame_map\": \"local_column_to_ensemble_reference_column\",\n"
     <<"  \"theta_zero\": \"forward\",\n  \"theta_pi\": \"backward\",\n  \"azimuth_periodic\": true,\n"
     <<"  \"input_polarization\": \"unpolarized\",\n  \"incident_total_intensity\": 1,\n"
     <<"  \"single_scattering_albedo\": 1,\n  \"stage\": \"incoherent\",\n"
     <<"  \"wavelength_nm\": "<<s.wavelength_nm<<",\n  \"interior_index\": "<<s.interior_index
     <<",\n  \"exterior_index\": "<<s.exterior_index<<",\n  \"query_frame_axis\": ";vector(f,s.grid.ki);
    f<<",\n  \"query_frame_e0\": ";vector(f,s.grid.frame.e0);
    f<<",\n  \"sampling_frame_columns\": [["<<s.grid.frame.e0.x<<','<<s.grid.frame.e1.x<<','<<s.grid.ki.x<<"],["
     <<s.grid.frame.e0.y<<','<<s.grid.frame.e1.y<<','<<s.grid.ki.y<<"],["
     <<s.grid.frame.e0.z<<','<<s.grid.frame.e1.z<<','<<s.grid.ki.z<<"]],\n"
     <<"  \"geometry\": {\"shape\": \"regular_hexagonal_prism\",\"body_long_axis\": [0,1,0],"
     <<"\"circumradius_mm\": "<<s.prism.circumradius_mm<<",\"full_length_mm\": "<<s.prism.length_mm
     <<",\"projected_area_mm2\": "<<area<<",\"projected_area_m2\": "<<area*1e-6
     <<",\"projected_area_method\": \"number_weighted_mean_of_convex_face_projections\"},\n"
     <<"  \"cross_sections\": {\"model\": \"geometric_nonabsorbing\",\"scattering_mm2\": "<<area
     <<",\"extinction_mm2\": "<<area<<",\"absorption_mm2\": 0,\"scattering_m2\": "<<area*1e-6
     <<",\"extinction_m2\": "<<area*1e-6<<",\"absorption_m2\": 0,\"wave_scattering_cross_section_computed\": false,"
     <<"\"assignment\": \"Csca=Cext=number_mean_Aproj; Cabs=0 (transport model)\"},\n"
     <<"  \"material\": {\"model\": \"caller_supplied_real_index_at_wavelength\","
     <<"\"wavelength_convention\": \"vacuum_nm\",\"absorption_model\": \"ignored_k_equals_zero\","
     <<"\"birefringence\": false,\"automatic_ice_dispersion\": false},\n"
     <<"  \"diffraction_policy\": {\"computed\": false,\"selected_output_includes_diffraction\": false},\n"
     <<"  \"inter_path_interference\": false,\n  \"tir_relative_phase_retained\": true,\n"
     <<"  \"orientation\": {\"representation\": \"weighted_discrete_SO3_number_measure\","
     <<"\"rotation\": \"active_Hamilton_quaternion_wxyz_body_to_ensemble\",\"rows\": "<<plan.nodes.size()
     <<",\"csv_sha256\": \""<<plan.source_sha256<<"\",\"input_weight_sum\": "<<plan.input_weight_sum
     <<",\"roundoff_weight_normalization\": "<<1/plan.input_weight_sum
     <<",\"plan_embedded\": false,\"retain_input_csv_for_reproduction\": true,\"symmetry_reduction\": false},\n"
     <<"  \"sampling\": {\"generator\": \"counter_splitmix64_open52_v1\",\"planned_rays\": "<<plan.total_samples
     <<",\"processed_rays\": "<<a.trace.rays<<",\"outgoing_samples\": "<<a.trace.outputs
     <<",\"batch_size\": "<<s.batch_size<<",\"max_output_records_per_batch\": "<<s.max_output_records
     <<",\"weight\": \"orientation_number_mass * Aproj / node_samples * path_power\"},\n"
     <<"  \"transport_audit\": {\"escaped_mm2\": "<<a.escaped_area.value<<",\"unresolved_mm2\": "<<a.unresolved_area.value
     <<",\"unresolved_fraction\": "<<a.unresolved_area.value/area<<",\"maximum_balance_error\": "<<a.trace.maximum_balance_error
     <<",\"maximum_interface_balance_error\": "<<a.trace.maximum_interface_balance_error
     <<",\"max_internal_hits\": "<<s.max_internal_hits<<",\"tail_tolerance_per_ray\": "<<s.tail_tolerance
     <<",\"max_unresolved_fraction\": "<<s.max_unresolved_fraction<<",\"max_balance_error\": "<<s.max_balance_error
     <<",\"cdf_normalized_to_resolved_escaped_mass\": true,\"status_counts\": {";
    for(unsigned i=0;i<8;++i){if(i)f<<',';quote(f,status_name(static_cast<TraceStatus>(i)));f<<':'<<a.trace.status_counts[i];}
    f<<"}},\n  \"histogram\": {\"theta_count\": "<<s.grid.nt<<",\"phi_count\": "<<s.grid.np
     <<",\"mass_mm2\": "<<a.grid_moments.mass<<",\"relative_mass_audit_error\": "<<a.histogram_relative_error
     <<",\"maximum_accepted_relative_mass_error\": "<<s.max_histogram_relative_error
     <<",\"g_rays\": "<<a.histogram.point_axial/a.histogram.point_mass
     <<",\"forward_axis_mass_mm2\": "<<a.histogram.forward_mass<<",\"backward_axis_mass_mm2\": "<<a.histogram.backward_mass
     <<",\"pole_policy\": \"axis_mass_uniform_in_azimuth_in_histogram_end_cap\",\"pole_roundoff_tolerance\": "<<pole_roundoff_tolerance
     <<",\"atomic_absorbed_addends\": "<<a.histogram.absorbed_addends<<",\"atomic_absorbed_addend_mass_mm2\": "<<a.histogram.absorbed_addend_mass
     <<",\"boundary_snapped_samples\": "<<a.histogram.boundary_snapped_samples
     <<",\"boundary_snap_rule\": \"64eps_coordinate_scaled_capped_at_1e-7_cell\""
     <<",\"atomic_diagnostic_is_error_bound\": false},\n"
     <<"  \"storage_processing\": {\"gaussian_sigma_degrees\": 0,\"mass_conservative_aggregation\": true,"
     <<"\"coarsening_tv\": "<<r.coarsening_tv<<",\"maximum_accepted_coarsening_tv\": "<<s.cdf_policy.max_coarsening_tv
     <<",\"aggregation_integral_relative_change\": "<<r.aggregation_relative_error
     <<",\"g_source\": "<<r.g_source<<",\"g_stored\": "<<r.g_stored<<"},\n"
     <<"  \"hg\": {\"g\": "<<r.g_stored<<",\"method\": \"first_moment_of_saved_cell_pdf\",\"target\": \"saved_cdf\","
     <<"\"cosine_convention\": \"dot(incident_propagation,outgoing_propagation)\"},\n"
     <<"  \"quality\": {\"source_coverage_certified\": false,\"angular_convergence_certified\": false,"
     <<"\"orientation_convergence_certified\": false,\"continuous_wave_solution_certified\": false,"
     <<"\"cdf_l1_mass_error\": "<<r.l1_error<<",\"cdf_maximum_cell_mass_error\": "<<r.maximum_cell_error
     <<",\"cdf_lost_positive_cells\": "<<r.lost_cells<<",\"cdf_lost_probability_mass\": "<<r.lost_mass
     <<",\"maximum_accepted_l1_mass_error\": "<<s.cdf_policy.max_l1_error
     <<",\"maximum_accepted_lost_probability_mass\": "<<s.cdf_policy.max_lost_mass<<",\"cdf_mass\": "<<r.cdf_mass<<"}\n}\n";
    f.close();
}
void write_failure(const std::filesystem::path& path,std::string_view reason,const CompositionAudit& a)
{
    auto f=output(path);f<<"{\n  \"complete\": false,\n  \"error\": ";quote(f,reason);
    f<<",\n  \"processed_rays\": "<<a.trace.rays<<",\n  \"orientations_done\": "<<a.orientations_done<<",\n  \"status_counts\": {";
    for(unsigned i=0;i<8;++i){if(i)f<<',';quote(f,status_name(static_cast<TraceStatus>(i)));f<<':'<<a.trace.status_counts[i];}
    f<<"}\n}\n";f.close();
}
} // namespace iceCrystal
