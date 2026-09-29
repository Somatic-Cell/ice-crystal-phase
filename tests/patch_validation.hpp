#pragma once

#include <rainbow/patch_geometry.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace rainbow::tests
{
inline void require_patch(const bool condition,const char* message)
{
    if(!condition) throw std::runtime_error(message);
}
inline bool close_patch(const float a,const float b)
{
    return std::isfinite(a) && std::isfinite(b)
        && std::abs(double(a)-double(b))<=1e-8+2e-5*std::max(std::abs(double(a)),std::abs(double(b)));
}

// Independent reference: integrate |P.(P_u x P_v)|/|P|^3 over [0,1]^2.
// This is composite midpoint quadrature in FP64, not the production atan2 formula.
// Used on smooth regular patches with a stated quadrature tolerance.
inline double integrated_patch_angle(const BilinearPatchGeometry& geometry,const unsigned count=160)
{
    struct D {double x,y,z;};
    const auto add=[](D a,D b){return D{a.x+b.x,a.y+b.y,a.z+b.z};};
    const auto sub=[](D a,D b){return D{a.x-b.x,a.y-b.y,a.z-b.z};};
    const auto mul=[](D a,double s){return D{a.x*s,a.y*s,a.z*s};};
    const auto dot=[](D a,D b){return a.x*b.x+a.y*b.y+a.z*b.z;};
    const auto cross=[](D a,D b){return D{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};};
    D p[4];
    for(unsigned i=0;i<4;++i) p[i]={geometry.corners[i].x,geometry.corners[i].y,geometry.corners[i].z};
    double integral=0;
    for(unsigned j=0;j<count;++j) for(unsigned i=0;i<count;++i)
    {
        const double u=(double(i)+0.5)/count,v=(double(j)+0.5)/count;
        const D q=add(add(mul(p[0],(1-u)*(1-v)),mul(p[1],u*(1-v))),
                      add(mul(p[2],(1-u)*v),mul(p[3],u*v)));
        const D du=add(mul(sub(p[1],p[0]),1-v),mul(sub(p[3],p[2]),v));
        const D dv=add(mul(sub(p[2],p[0]),1-u),mul(sub(p[3],p[1]),u));
        const double r2=dot(q,q);
        integral+=std::abs(dot(q,cross(du,dv)))/(r2*std::sqrt(r2));
    }
    return integral/(double(count)*count);
}

inline BilinearPatchGeometry square_patch()
{
    return {{Vec3{-1,-1,1}.normalized(),Vec3{1,-1,1}.normalized(),
             Vec3{-1,1,1}.normalized(),Vec3{1,1,1}.normalized()}};
}

inline std::vector<OutgoingVertex> synthetic_patch_vertices()
{
    // Four 2x2 families: regular, folded, missing corner, no outgoing corners.
    std::vector<OutgoingVertex> vertices(16);
    const auto geometry=square_patch();
    for(unsigned p=0;p<4;++p) for(unsigned k=0;k<4;++k)
    {
        auto& v=vertices[p*4+k];
        v.status=VertexStatus::Valid;
        v.direction_drop=geometry.corners[k];
        v.basis_x={1,0,0};
        v.field={{1,0},{0,0}};
    }
    std::swap(vertices[5].direction_drop,vertices[7].direction_drop);
    vertices[8].status=VertexStatus::Miss;
    for(unsigned k=12;k<16;++k) vertices[k].status=VertexStatus::TotalInternalReflection;
    return vertices;
}

inline void verify_patch_buffer(
    std::span<const OutgoingVertex> vertices,const PatchBuildLayout layout,
    std::span<const PatchCellStatus> statuses,std::span<const OutgoingPatch> patches,
    std::span<const PatchAabb> aabbs)
{
    require_patch(statuses.size()==layout.cell_count,"Cell status count mismatch.");
    require_patch(aabbs.size()==patches.size(),"AABB count mismatch.");
    std::size_t compact_index=0;
    for(std::uint32_t id=0;id<layout.cell_count;++id)
    {
        const auto reference=PatchConstruction::make(vertices.data(),layout,id);
        require_patch(statuses[id]==reference.patch.status,"CPU/GPU patch classification mismatch.");
        if(!has_patch_geometry(statuses[id])) continue;
        require_patch(compact_index<patches.size(),"Missing compacted patch.");
        const auto& p=patches[compact_index];
        const auto& b=aabbs[compact_index];
        require_patch(p.patch_id==id && p.status==statuses[id],"Stable patch ID/order mismatch.");
        for(unsigned k=0;k<4;++k)
            require_patch(p.vertex_indices[k]==reference.patch.vertex_indices[k],"Corner mapping mismatch.");
        require_patch(p.incident_area_drop2==layout.incident_area_drop2,"Incident cell area changed.");
        if(p.has_regular_spherical_map())
            require_patch(close_patch(p.signed_solid_angle_sr,reference.patch.signed_solid_angle_sr),"Solid angle mismatch.");
        else require_patch(std::isnan(p.signed_solid_angle_sr),"Nonregular patch received a fabricated solid angle.");
        const auto& r=reference.aabb;
        require_patch(b.min_x==r.min_x && b.min_y==r.min_y && b.min_z==r.min_z
                   && b.max_x==r.max_x && b.max_y==r.max_y && b.max_z==r.max_z,"AABB mismatch.");
        BilinearPatchGeometry geometry{};
        for(unsigned k=0;k<4;++k) geometry.corners[k]=vertices[p.vertex_indices[k]].direction_drop;
        for(unsigned v=0;v<=4;++v) for(unsigned u=0;u<=4;++u)
        {
            const auto point=geometry.evaluate(float(u)*0.25f,float(v)*0.25f);
            require_patch(point.x>=b.min_x && point.x<=b.max_x
                       && point.y>=b.min_y && point.y<=b.max_y
                       && point.z>=b.min_z && point.z<=b.max_z,"Bilinear surface escaped its AABB.");
        }
        ++compact_index;
    }
    require_patch(compact_index==patches.size(),"Extra compacted patches.");
}
} // namespace rainbow::tests
