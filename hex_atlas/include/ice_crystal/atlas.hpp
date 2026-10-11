#pragma once
#include <ice_crystal/hex_trace.hpp>
#include <ice_crystal/phase_math.hpp>
#include <cstddef>
#include <filesystem>
#include <span>
#include <vector>

namespace iceCrystal::atlas
{
struct AreaRay
{
    EntrySample entry{};
    double projected_weight_mm2=0;
};
// Tensor Gauss-Legendre (2 x 3) in each equal-area parameter cell. No RNG.
std::vector<AreaRay> make_area_rule(const HexPrism&,Vec3 ki,std::uint32_t subdivisions);
struct Settings
{
    HexPrism prism{};
    double interior_index=0,exterior_index=1,wavelength_nm=0;
    std::uint32_t incident_theta_intervals=0,incident_phi_intervals=0;
    std::uint32_t theta_count=0,phi_count=0,area_subdivisions=0,max_internal_hits=2048;
    std::size_t batch_rays=64,max_work_bytes=512ull*1024*1024;
    double tail_tolerance=1e-12,max_unresolved_fraction=1e-10,max_balance_error=1e-10;
};
struct Node
{
    double theta=0,phi=0;
    TraceSettings trace{};
    std::vector<AreaRay> area_rays;
    std::vector<double> u_edges;
    PhaseGridView grid{};
    // Rebind grid.u_edges after a move/copy; consumers use view(), not grid.u_edges.
    PhaseGridView view() const noexcept {auto v=grid;v.u_edges=u_edges.data();return v;}
};
Node make_node(const Settings&,std::uint32_t theta_index,std::uint32_t phi_index);
void validate(const Settings&);
void validate_node(const Settings&,const Node&);
struct Result
{
    std::vector<double> mass_mm2;
    TraceSummary audit{};
    double projected_area_mm2=0,escaped_mm2=0,unresolved_mm2=0;
    double histogram_relative_error=0;
    Sum escaped_accumulator{},unresolved_accumulator{};
};
Result trace_cpu(const Settings&,const Node&);
// CPU and CUDA call the identical audit code on per-ray accounting.
void add_audits(Result&,std::span<const RayAudit>,std::span<const AreaRay>,double tolerance);
void finish(Result&,const Settings&,const Node&);

struct MassSink
{
    PhaseGridView grid{};
    Sum* mass=nullptr; // bin-major, ray-local columns
    std::uint64_t stride=1,column=0;
    double weight=0;
    ICE_HD bool add(std::uint64_t bin,double w) noexcept
    {
        if(!(w>0)||!finite_value(w))return false;
        mass[bin*stride+column].add(w);
        return finite_value(mass[bin*stride+column].value);
    }
    ICE_HD bool push(const OutgoingSample& s) noexcept
    {
        const auto b=phase_bin(grid,s.direction);
        const double w=weight*s.power_fraction;
        if(!b.valid||!(w>0)||!finite_value(w))return false;
        if(!b.pole)return add(b.index,w);
        const auto row=b.pole>0?0:grid.nt-1;
        for(std::uint32_t j=0;j<grid.np;++j)
            if(!add(std::uint64_t(row)*grid.np+j,w/grid.np))return false;
        return true;
    }
};
struct BatchParams
{
    TraceSettings trace{};
    PhaseGridView grid{};
    const AreaRay* rays=nullptr;
    Sum* mass=nullptr;
    RayAudit* audits=nullptr;
    std::uint64_t first=0,count=0;
};
static_assert(std::is_trivially_copyable_v<AreaRay>);
static_assert(std::is_trivially_copyable_v<BatchParams>);
}
