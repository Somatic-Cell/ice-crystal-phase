#include <optix.h>
#include <cstdio>
#include <rainbow/raindrop_paths.hpp>

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
        // 前回と同じ位置：Unresolved へ変更する前の結果を表示する．
        if((query.hit.diagnostics & rainbow::IntersectionFallback) != 0u
           || query.hit.code == rainbow::BoundaryCode::Unresolved
           || (query.hit.code == rainbow::BoundaryCode::Hit && !query.accepted))
        {
            printf("[drop/query] vertex=%u inside=%u code=%u closest_hit=%u "
                   "diagnostics=%u distance=%.9g\n",
                optixGetLaunchIndex().x,
                static_cast<unsigned>(ray.starts_inside),
                static_cast<unsigned>(query.hit.code),
                static_cast<unsigned>(query.accepted),query.hit.diagnostics,
                static_cast<double>(query.hit.distance));
        }
        if(query.hit.code==rainbow::BoundaryCode::Hit && !query.accepted)
            query.hit.code=rainbow::BoundaryCode::Unresolved;
        return query.hit;
    }
};
}

extern "C" __global__ void __raygen__raindrop_trace()
{
    const unsigned index=optixGetLaunchIndex().x;
    if(index>=raindrop_trace_params.config->vertex_count()) return;
    OptixBoundaryIntersector intersector;
    rainbow::RaindropPathTracer::trace_vertex(*raindrop_trace_params.config,index,
        intersector,raindrop_trace_params.vertices);
}

extern "C" __global__ void __intersection__raindrop()
{
    auto* query=query_from_payload();
    rainbow::RaindropIntersector intersector{raindrop_trace_params.config->shape};
    // 原因特定用．曲面・許容誤差・反復上限は元の実装のまま．
    rainbow::detail::BoundarySolveDiagnostics diagnostic{};
    query->hit=intersector.intersect(query->ray,&diagnostic);
    if(query->hit.code==rainbow::BoundaryCode::Unresolved)
    {
        const unsigned index=optixGetLaunchIndex().x;
        printf("[drop/root-stage] vertex=%u inside=%u fp32_code=%u fp32_reason=%u "
               "used_fp64=%u fp64_code=%u fp64_reason=%u\n",
            index,static_cast<unsigned>(query->ray.starts_inside),
            static_cast<unsigned>(diagnostic.fp32.result),
            static_cast<unsigned>(diagnostic.fp32.reason),
            static_cast<unsigned>(diagnostic.used_fp64),
            static_cast<unsigned>(diagnostic.fp64.result),
            static_cast<unsigned>(diagnostic.fp64.reason));

        if(diagnostic.used_fp64)
        {
            const auto& d=diagnostic.fp64;
            printf("[drop/root-stop] vertex=%u reason=%u visits=%u iterations=%u "
                   "pending=%u lo=%.17g hi=%.17g width=%.17g tol=%.17g\n",
                index,static_cast<unsigned>(d.reason),d.visits,d.iterations,
                d.pending_segments,d.lower,d.upper,d.upper-d.lower,d.tolerance);
            printf("[drop/root-values] vertex=%u has_f=%u f_lo=%.17g f_hi=%.17g "
                   "has_bounds=%u F=[%.17g,%.17g] D=[%.17g,%.17g]\n",
                index,static_cast<unsigned>(d.has_endpoint_values),d.f_lower,d.f_upper,
                static_cast<unsigned>(d.has_bounds),d.bounds.value.lo,d.bounds.value.hi,
                d.bounds.derivative.lo,d.bounds.derivative.hi);
        }

        // 既知の失敗レイを一つだけ詳しく調べる．他のレイの処理は変えない．
        if(index==3884u && !query->ray.starts_inside)
        {
            const auto o=query->ray.physical_origin;
            const auto w=query->ray.direction;
            printf("[drop/root-ray] vertex=%u o=(%.9g,%.9g,%.9g) w=(%.9g,%.9g,%.9g)\n",
                index,static_cast<double>(o.x),static_cast<double>(o.y),static_cast<double>(o.z),
                static_cast<double>(w.x),static_cast<double>(w.y),static_cast<double>(w.z));
            // 入力レイに由来する値なので，定数だけのコンパイル時評価ではない．
            const double probe=static_cast<double>(w.x);
            const double down=rainbow::detail::next_down(probe);
            const double up=rainbow::detail::next_up(probe);
            printf("[drop/fp64-step] vertex=%u x=%.17g down=%.17g up=%.17g "
                   "down_step=%.17g up_step=%.17g\n",
                index,probe,down,up,probe-down,up-probe);
        }
    }
    if(query->hit.code==rainbow::BoundaryCode::Hit)
    {
        const float t=query->hit.distance-query->ray.start_distance;
        // closest-hit が採用した場合だけ accepted=true にする．
        const float ray_t_min=optixGetRayTmin();
        const float ray_t_max=optixGetRayTmax();
        const bool was_reported=optixReportIntersection(t,0);
        if(!was_reported)
        {
            printf("[drop/report_rejected] vertex=%u inside=%u t=%.9g t_min=%.9g t_max=%.9g\n",
                optixGetLaunchIndex().x,static_cast<unsigned>(query->ray.starts_inside),
                static_cast<double>(t),static_cast<double>(ray_t_min),static_cast<double>(ray_t_max));
        }
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
