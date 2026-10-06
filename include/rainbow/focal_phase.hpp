#pragma once
#include <rainbow/wave_optics_data.hpp>
#include <cmath>

namespace rainbow
{
struct FocalPhaseEstimate
{
    double radial_derivative = 0;
    std::uint32_t family = 0, extra_quarter_turn = 0, quarter_turns = 0;
};
// Sec.4.1 "Focal Lines". b=U^2+V^2 uses EMITTING-plane coordinates.
// Discretization: radial derivative of bilinear corner theta at cell center.
// Regular: cell center, as before. Folded extension: one representative per
// connected branch. In both cases phase is independent of query-hit position.
struct FocalPhase
{
    [[nodiscard]] HOST_DEVICE static bool layout_valid(const FocalPhaseConfig& c) noexcept
    {
        if(c.grid_width < 2 || c.grid_height < 2 || !(c.grid_half_extent > 0)
           || !(c.grid_half_extent <= 0x1.fffffep127f)) return false;
        if(std::uint64_t(c.grid_width)*c.grid_height > 0xffffffffull/4u
           || std::uint64_t(c.grid_width-1u)*(c.grid_height-1u) > 0xffffffffull/4u) return false;
        for(unsigned j=0; j<4; ++j) if(c.quarter_turn_offsets[j] > 3) return false;
        return true;
    }
    [[nodiscard]] HOST_DEVICE static std::uint32_t estimate(
        const OutgoingVertex* vertices, std::uint32_t vertex_count,
        const OutgoingPatch& patch, const Vec3 incident, const FocalPhaseConfig& c,
        FocalPhaseEstimate& output) noexcept
    {
        return estimate_impl(vertices, vertex_count, patch, incident, c, 0.5, 0.5, true, output);
    }
    // Branch representative in the ORIGINAL incident cell. This is not J.
    [[nodiscard]] HOST_DEVICE static std::uint32_t estimate_at(
        const OutgoingVertex* vertices, std::uint32_t vertex_count,
        const OutgoingPatch& patch, const Vec3 incident, const FocalPhaseConfig& c,
        double u, double v, FocalPhaseEstimate& output) noexcept
    {
        return estimate_impl(vertices, vertex_count, patch, incident, c, u, v, false, output);
    }
    [[nodiscard]] HOST_DEVICE static OpticalField64 apply(OpticalField64 a,std::uint32_t q) noexcept
    {
        // exp(+i*pi*q/2), exact signs/swaps, no OPL or trigonometric re-evaluation.
        switch(q&3u)
        {
        case 0:return a;
        case 1:return {-a.s_imag,a.s_real,-a.p_imag,a.p_real};
        case 2:return {-a.s_real,-a.s_imag,-a.p_real,-a.p_imag};
        default:return {a.s_imag,-a.s_real,a.p_imag,-a.p_real};
        }
    }
private:
    [[nodiscard]] HOST_DEVICE static std::uint32_t estimate_impl(
        const OutgoingVertex* vertices, std::uint32_t vertex_count,
        const OutgoingPatch& patch, const Vec3 incident, const FocalPhaseConfig& c,
        double u, double v, bool at_center, FocalPhaseEstimate& output) noexcept
    {
        if(!(u>=0&&u<=1&&v>=0&&v<=1) || !vertices || !layout_valid(c) || !incident.is_finite()) return FocalInvalidGeometry;
        const auto n = c.grid_width*c.grid_height;
        const auto nx = c.grid_width-1u;
        const auto cells = nx*(c.grid_height-1u);
        if(patch.patch_id >= 4u*cells || vertex_count < 4u*n) return FocalInvalidGeometry;
        const auto family=patch.patch_id/cells, cell=patch.patch_id%cells;
        const auto ix=cell%nx, iy=cell/nx;
        const auto first=family*n+iy*c.grid_width+ix;
        const std::uint32_t expected[4]={first,first+1u,first+c.grid_width,first+c.grid_width+1u};
        const auto wi=incident.cast<double>();
        if(!(wi.dot(wi)>0)) return FocalInvalidGeometry;
        double theta[4]{};
        for(unsigned j=0;j<4;++j)
        {
            if(patch.vertex_indices[j]!=expected[j]) return FocalInvalidGeometry;
            const auto& a=vertices[expected[j]];
            if(a.status!=VertexStatus::Valid || !a.direction_drop.is_finite()) return FocalInvalidGeometry;
            const auto wo=a.direction_drop.cast<double>();
            if(!(wo.dot(wo)>0)) return FocalInvalidGeometry;
            theta[j]=::atan2(wi.cross(wo).length(),wi.dot(wo));
            if(!finite(theta[j])) return FocalInvalidGeometry;
        }
        const double e=double(c.grid_half_extent);
        const double du=2*e/double(c.grid_width-1u), dv=2*e/double(c.grid_height-1u);
        const double U=(2*(double(ix)+u)/double(c.grid_width-1u)-1)*e;
        const double V=(2*(double(iy)+v)/double(c.grid_height-1u)-1)*e;
        const double r2=::fma(U,U,V*V);
        if(!(r2>0)) return FocalDerivativePending;
        const double gu=at_center ? ((theta[1]-theta[0])+(theta[3]-theta[2]))/(2*du)
                                  : ((1-v)*(theta[1]-theta[0])+v*(theta[3]-theta[2]))/du;
        const double gv=at_center ? ((theta[2]-theta[0])+(theta[3]-theta[1]))/(2*dv)
                                  : ((1-u)*(theta[2]-theta[0])+u*(theta[3]-theta[1]))/dv;
        const double numerator=::fma(U,gu,V*gv);
        // Arithmetic cancellation guard, NOT a bound on input/geometric error.
        const double guard=0x1p-46*(::fabs(U*gu)+::fabs(V*gv)+1);
        if(!finite(numerator) || !finite(r2)) return FocalInvalidGeometry;
        if(::fabs(numerator)<=guard) return FocalDerivativePending;
        FocalPhaseEstimate a{};
        a.family=family; a.radial_derivative=numerator/(2*r2);
        a.extra_quarter_turn=numerator>0 ? 1u : 0u;
        a.quarter_turns=(c.quarter_turn_offsets[family]+a.extra_quarter_turn)&3u;
        if(!finite(a.radial_derivative)) return FocalInvalidGeometry;
        output=a; return FocalNone;
    }
    [[nodiscard]] HOST_DEVICE static bool finite(double x) noexcept
    {return x>=-0x1.fffffffffffffp1023 && x<=0x1.fffffffffffffp1023;}
};
} // namespace rainbow
