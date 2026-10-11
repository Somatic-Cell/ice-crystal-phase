#include <ice_crystal/ice_optics.hpp>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <stdexcept>

namespace iceCrystal::optics
{
namespace
{
struct Entry { double nm, n, k; };
#include "ice_optics_table.inc"
}
IceState evaluate_ice_wb2008(double wavelength_nm)
{
    if(!std::isfinite(wavelength_nm) || wavelength_nm < 200 || wavelength_nm > 1100)
        throw std::invalid_argument("Warren-Brandt 2008 ice: supported vacuum wavelengths are 200..1100 nm; no extrapolation.");
    const auto* b=std::lower_bound(std::begin(ice_table),std::end(ice_table),wavelength_nm,
        [](const Entry& e,double x){return e.nm<x;});
    if(b==std::end(ice_table))throw std::logic_error("Missing ice optical bracket.");
    if(b->nm==wavelength_nm)return {wavelength_nm,b->n,b->k,b->nm,b->nm,0};
    const auto& a=*(b-1);
    const double t=std::log(wavelength_nm/a.nm)/std::log(b->nm/a.nm);
    // The author specifies n versus log(lambda), and log(k) versus log(lambda).
    const double n=a.n+t*(b->n-a.n);
    const double k=std::exp(std::log(a.k)+t*(std::log(b->k)-std::log(a.k)));
    if(!(n>0 && k>0 && std::isfinite(n) && std::isfinite(k)))
        throw std::runtime_error("Invalid ice optical state.");
    return {wavelength_nm,n,k,a.nm,b->nm,t};
}
void write_material_json(std::ostream& o,const IceState& s,double exterior,bool ice)
{
    o<<"{\"model\":\""<<(ice?"ice_ih_warren_brandt_2008_real_only":"constant")<<"\",";
    o<<"\"wavelength_convention\":\"vacuum_nm\",\"wavelength_nm\":"<<s.wavelength_nm
     <<",\"wavelength_requested_nm\":"<<s.wavelength_nm<<",\"wavelength_evaluated_nm\":"<<s.wavelength_nm
     <<",\"interior_index\":"<<s.n<<",\"exterior_index\":"<<exterior
     <<",\"interior_index_evaluated_fp64\":"<<s.n<<",\"interior_index_solver_fp64\":"<<s.n
     <<",\"automatic_ice_dispersion\":"<<(ice?"true":"false")
     <<",\"absorption\":false,\"birefringence\":false,\"absorption_model\":\"ignored_k_equals_zero\",\"kappa_solver\":0";
    if(ice)
        o<<",\"source_doi\":\"10.1029/2007JD009744\",\"reference_state\":\"ice_Ih_near_melting_compilation\""
         <<",\"source_url\":\"https://www.atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat\""
         <<",\"table_subset_sha256\":\"71235dd08b3f65fc5fdf638101577d4da12af9cd725da21fa1677c0f91bfddcb\",\"supported_nm\":[200,1100]"
         <<",\"n_interpolation\":\"linear_n_in_log_vacuum_wavelength\",\"kappa_interpolation\":\"linear_log_kappa_in_log_vacuum_wavelength\""
         <<",\"temperature_dependence_computed\":false,\"kappa_evaluated_not_used\":"<<s.kappa
         <<",\"bracket_nm\":["<<s.lower_nm<<','<<s.upper_nm<<"],\"log_wavelength_fraction\":"<<s.log_wavelength_fraction;
    o<<'}';
}
} // namespace iceCrystal::optics
