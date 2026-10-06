#include <optix.h>
#include <rainbow/raindrop_paths.hpp>

// GPU 側から見えるパラメータ
extern "C"
{
__constant__ rainbow::RaindropTraceParams raindrop_trace_params;
}

namespace
{
struct BoundaryQuery
{
    rainbow::BoundaryRay ray;
    rainbow::BoundaryHit hit;
    bool accepted=false;
};

__forceinline__ __device__ BoundaryQuery* query_from_payload()
{
    const unsigned long long address=static_cast<unsigned long long>(optixGetPayload_0())
        |(static_cast<unsigned long long>(optixGetPayload_1())<<32);
    return reinterpret_cast<BoundaryQuery*>(address);
}

struct OptixBoundaryIntersector
{
    __device__ rainbow::BoundaryHit intersect(const rainbow::BoundaryRay& ray) const noexcept
    {
        BoundaryQuery query{ray,{},false};
        const auto address=reinterpret_cast<unsigned long long>(&query);
        unsigned lo=static_cast<unsigned>(address),hi=static_cast<unsigned>(address>>32);
        const rainbow::Vec3 origin=ray.physical_origin.at(ray.direction,ray.start_distance);
        const float far=8.0f*raindrop_trace_params.config->shape.outer_radius+4.0f;
        // Offset されたレイは GAS の traversal 専用．交差計算は query 内の物理 ray を使う．
        // 一つの保守的 AABB が曲面全体を包み，数値 origin の丸めで候補を失わない．
        optixTrace(static_cast<OptixTraversableHandle>(raindrop_trace_params.traversable),
            make_float3(origin.x,origin.y,origin.z),
            make_float3(ray.direction.x,ray.direction.y,ray.direction.z),
            0.0f,far,0.0f,255,OPTIX_RAY_FLAG_DISABLE_ANYHIT,
            0,1,0,lo,hi);
        if(query.hit.code==rainbow::BoundaryCode::Hit && !query.accepted)
            query.hit.code=rainbow::BoundaryCode::Unresolved;
        return query.hit;
    }
};
}

extern "C" __global__ void __raygen__raindrop_trace_v2()
{
    const unsigned index=optixGetLaunchIndex().x;
    if(index>=raindrop_trace_params.config->vertex_count()) return;
    OptixBoundaryIntersector intersector;
    if(raindrop_trace_params.second_input_fields)
        rainbow::RaindropPathTracer::trace_vertex_unpolarized(*raindrop_trace_params.config,index,
            intersector,raindrop_trace_params.vertices,raindrop_trace_params.second_input_fields);
    else
        rainbow::RaindropPathTracer::trace_vertex(*raindrop_trace_params.config,index,
            intersector,raindrop_trace_params.vertices);
}

extern "C" __global__ void __intersection__raindrop()
{
    auto* query=query_from_payload();
    rainbow::RaindropIntersector intersector{raindrop_trace_params.config->shape};
    // Normal execution passes no RootSolverDiagnostics. The corrected typed
    // nextafter calls and all solver tolerances remain in the existing header.
    query->hit=intersector.intersect(query->ray);
    if(query->hit.code==rainbow::BoundaryCode::Hit)
    {
        const float t=query->hit.distance-query->ray.start_distance;
        // closest-hit が採用した場合だけ accepted=true にする．
        static_cast<void>(optixReportIntersection(t,0));
    }
}
extern "C" __global__ void __closesthit__raindrop()
{
    query_from_payload()->accepted=true;
}
extern "C" __global__ void __miss__raindrop()
{
    // IS が Unresolved を書いた場合に，単なる miss で上書きしない．
}
