#pragma once

#include <rainbow/raindrop_trace_data.hpp>
#include <limits>
#include <stdexcept>

namespace rainbow
{
// 実行条件．I/O / AS / device pointer を持たないホスト側設定．
struct RaindropSettings
{
    float radius_mm=0.4f;
    float wavelength_nm=700.0f;
    float interior_index=1.3314f; // paper p.4 の 700nm の例．他波長では利用者が値を与える．
    Vec3 incident_direction{1,0,0};
    Field32 incident_field{{1,0},{0,0}};
    std::uint32_t grid_width=129,grid_height=129;
    bool force_sphere=false;

    [[nodiscard]] RaindropTraceConfig make_config() const
    {
        RaindropTraceConfig c{};
        if(!RaindropShape::try_make(radius_mm,force_sphere,c.shape))
            throw std::invalid_argument("radius_mm must be in (0,3]; no coefficient extrapolation is implemented.");
        if(!(wavelength_nm>=380.0f && wavelength_nm<=830.0f))
            throw std::invalid_argument("wavelength_nm must be in [380,830] for this implementation.");
        if(!(interior_index>=1.0f && interior_index<=2.0f))
            throw std::invalid_argument("interior_index must be in [1,2]; exterior_index is 1.");
        if(grid_width<2 || grid_height<2 || std::uint64_t(grid_width)*grid_height>std::uint64_t((std::numeric_limits<std::uint32_t>::max)())/4)
            throw std::invalid_argument("Invalid grid dimensions or 32-bit path-major index overflow.");
        if(!incident_direction.is_finite() || !std::isfinite(incident_direction.length()) || !(incident_direction.length()>0))
            throw std::invalid_argument("incident_direction must be a finite nonzero vector.");
        if(!std::isfinite(incident_field.x.real)||!std::isfinite(incident_field.x.imag)
           ||!std::isfinite(incident_field.y.real)||!std::isfinite(incident_field.y.imag))
            throw std::invalid_argument("incident_field must be finite.");
        c.radius_mm=radius_mm; c.wavelength_nm=wavelength_nm; c.interior_index=interior_index;
        c.incident_direction=incident_direction.normalized();
        c.incident_basis_x=TransverseFrame::from_direction(c.incident_direction).e0;
        c.incident_field=incident_field;
        c.grid_width=grid_width; c.grid_height=grid_height;
        c.grid_half_extent=c.shape.outer_radius*1.01f;
        c.reference_distance=c.outgoing_reference_distance=2.0f*c.shape.outer_radius;
        return c;
    }
};
}
