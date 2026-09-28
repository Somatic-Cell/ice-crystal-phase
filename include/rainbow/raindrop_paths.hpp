#pragma once

#include <rainbow/raindrop_trace_data.hpp>
#include <rainbow/dielectric_interface.hpp>
#include <rainbow/optical_path.hpp>

namespace rainbow
{

// 一つの格子 vertex から R を放出し，T を継続．各内部境界で T を放出して R を継続する．
// CPU/OptiX は境界問い合わせだけを差し替える．アルゴリズム・phase/field 演算は共有する．
struct RaindropPathTracer
{
    template<class Intersector>
    HOST_DEVICE static void trace_vertex(
        const RaindropTraceConfig& config,const std::uint32_t index,
        Intersector& intersector,OutgoingVertex* vertices) noexcept
    {
        const std::uint32_t count=config.vertex_count();
        for(unsigned family=0;family<4;++family) vertices[family*count+index]=OutgoingVertex{};
        const std::uint32_t ix=index%config.grid_width,iy=index/config.grid_width;
        const float u=(2.0f*static_cast<float>(ix)/static_cast<float>(config.grid_width-1)-1.0f)*config.grid_half_extent;
        const float v=(2.0f*static_cast<float>(iy)/static_cast<float>(config.grid_height-1)-1.0f)*config.grid_half_extent;
        const Vec3 w=config.incident_direction;
        const TransverseFrame frame{config.incident_basis_x,w.cross(config.incident_basis_x).normalized()};
        PolarizedRay active{w,frame,config.incident_field};
        const Vec3 grid_point=(frame.e0*u+frame.e1*v)-w*config.reference_distance;
        const BoundaryHit entry=intersector.intersect({grid_point,w,0.0f,false});
        std::uint32_t diagnostics=entry.diagnostics;
        if(entry.code!=BoundaryCode::Hit)
        {
            mark_remaining(vertices,count,index,0,
                entry.code==BoundaryCode::Miss?VertexStatus::Miss:VertexStatus::UnresolvedIntersection,diagnostics);
            return;
        }
        OpticalPathAccumulator path=OpticalPathAccumulator::make(config.radius_mm,config.wavelength_nm);
        const FloatPair before=FloatPair{config.reference_distance,0}+FloatPair::dot(entry.position,w);
        if(!path.add_distance(before,config.exterior_index))
        {mark_remaining(vertices,count,index,0,VertexStatus::PhaseOverflow,diagnostics);return;}
        const auto branches=DielectricInterface::split(active,entry.outward_normal,config.exterior_index,config.interior_index);
        diagnostics|=branches.diagnostics;
        if(!branches.is_valid)
        {mark_remaining(vertices,count,index,0,VertexStatus::InvalidInterface,diagnostics);return;}
        emit(vertices[index],branches.reflected,entry.position,path,config,diagnostics);
        if(!branches.has_transmission)
        {mark_remaining(vertices,count,index,1,VertexStatus::TotalInternalReflection,diagnostics);return;}
        active=branches.transmitted;
        Vec3 previous_position=entry.position;
        for(unsigned family=1;family<4;++family)
        {
            // along-ray self-hit exclusion．query origin と physical point を混同しない．
            // 全光路は previous_position/実交点から計算し，この offset を足さない．
            constexpr float surface_offset=0x1p-17f;
            const BoundaryHit hit=intersector.intersect({previous_position,active.direction,surface_offset,true});
            diagnostics|=hit.diagnostics;
            if(hit.code!=BoundaryCode::Hit)
            {mark_remaining(vertices,count,index,family,VertexStatus::UnresolvedIntersection,diagnostics);return;}
            if(!path.add_segment(previous_position,hit.position,config.interior_index))
            {mark_remaining(vertices,count,index,family,VertexStatus::PhaseOverflow,diagnostics);return;}
            const auto split=DielectricInterface::split(active,-hit.outward_normal,config.interior_index,config.exterior_index);
            diagnostics|=split.diagnostics;
            if(!split.is_valid)
            {mark_remaining(vertices,count,index,family,VertexStatus::InvalidInterface,diagnostics);return;}
            if(split.has_transmission)
                emit(vertices[family*count+index],split.transmitted,hit.position,path,config,diagnostics);
            else
            {
                auto& vertex=vertices[family*count+index];
                vertex.status=VertexStatus::TotalInternalReflection;
                vertex.diagnostics=diagnostics;
            }
            // Brewster/小振幅/零電場でも枝刈りしない．TIR でも内部反射を継続する．
            active=split.reflected;
            previous_position=hit.position;
        }
    }
private:
    HOST_DEVICE static void mark_remaining(OutgoingVertex* vertices,const std::uint32_t count,
        const std::uint32_t index,const unsigned first,const VertexStatus status,const std::uint32_t diagnostics) noexcept
    {
        for(unsigned p=first;p<4;++p)
        {vertices[p*count+index].status=status;vertices[p*count+index].diagnostics=diagnostics;}
    }
    HOST_DEVICE static void emit(OutgoingVertex& vertex,const PolarizedRay& outgoing,
        const Vec3 physical_point,OpticalPathAccumulator path,const RaindropTraceConfig& config,
        const std::uint32_t diagnostics) noexcept
    {
        const FloatPair after=FloatPair{config.outgoing_reference_distance,0}
            -FloatPair::dot(physical_point,outgoing.direction);
        vertex.diagnostics=diagnostics;
        const auto& f=outgoing.field;
        if(!(::fabsf(f.x.real)<=0x1.fffffep127f && ::fabsf(f.x.imag)<=0x1.fffffep127f
             && ::fabsf(f.y.real)<=0x1.fffffep127f && ::fabsf(f.y.imag)<=0x1.fffffep127f))
        {vertex.status=VertexStatus::InvalidInterface;return;}
        if(!path.add_distance(after,config.exterior_index))
        {vertex.status=VertexStatus::PhaseOverflow;return;}
        vertex.position_drop=physical_point;
        vertex.direction_drop=outgoing.direction;
        vertex.basis_x=outgoing.frame.e0;
        vertex.field=outgoing.field;
        vertex.optical_cycles=path.phase;
        vertex.status=VertexStatus::Valid;
    }
};
}
