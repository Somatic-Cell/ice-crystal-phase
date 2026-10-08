#pragma once
#include <ice_crystal/phase_math.hpp>
#include <filesystem>
#include <span>
#include <vector>

namespace iceCrystal
{
struct PhaseGrid
{
    std::uint32_t nt=0,np=0;
    Vec3 ki{};
    Frame frame{};
    std::vector<double> edges;
    [[nodiscard]] PhaseGridView view() const noexcept {return {edges.data(),ki,frame,nt,np};}
};
PhaseGrid make_phase_grid(std::uint32_t nt,std::uint32_t np,Vec3 ki);
void validate_phase_grid(const PhaseGrid&);

// Nonnegative integrated cross-section contributions, never PDF samples.
class PhaseMassAccumulator final
{
public:
    explicit PhaseMassAccumulator(PhaseGrid);
    void add(Vec3 outgoing_reference,double weight_mm2);
    void add(std::span<const OutgoingSample>,const Rotation&,double area_weight_per_ray);
    [[nodiscard]] std::vector<double> masses() const;
    [[nodiscard]] HistogramStats statistics() const;
    [[nodiscard]] const PhaseGrid& grid() const noexcept {return grid_;}
private:
    PhaseGrid grid_;
    std::vector<Sum> mass_;
    Sum point_mass_,point_axial_,forward_,backward_;
    std::uint64_t count_=0,snapped_=0;
};
HistogramMoments measure_mass(const PhaseGrid&,std::span<const double>);
std::vector<double> mass_to_density(const PhaseGrid&,std::span<const double>);

struct CdfPolicy
{
    double max_l1_error=1e-10,max_lost_mass=1e-12,max_coarsening_tv=1;
};
struct CdfReport
{
    double g_source=0,g_stored=0,coarsening_tv=0,cdf_mass=0;
    double l1_error=0,lost_mass=0,maximum_cell_error=0,aggregation_relative_error=0;
    std::uint64_t lost_cells=0;
};
struct HostCdf
{
    PhaseGrid grid;
    std::vector<double> phi,theta;
    CdfReport report;
    void write_arrays(const std::filesystem::path&) const;
};
// CPU reference: conservative mass aggregation + ordered binary64 scans.
// CUDA production path instead calls the EXISTING rainbow::PhaseCdf.
HostCdf build_host_cdf(const PhaseGrid&,std::span<const double> masses,
    std::uint32_t out_nt,std::uint32_t out_np,const CdfPolicy& policy={});
void validate_cdf_policy(const CdfPolicy&);
} // namespace iceCrystal
