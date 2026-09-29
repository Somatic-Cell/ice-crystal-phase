#include <optix.h>
#include <rainbow/patch_query_launch_params.hpp>

extern "C"
{
__constant__ rainbow::PatchQueryLaunchParams patch_query_params;
}

namespace
{
struct QueryState
{
    rainbow::PatchHitCollector collector{};
    rainbow::BilinearRayHits pending{};
    std::uint32_t compact_index=0;
    bool consumed=false;
};

__forceinline__ __device__ QueryState* query_state()
{
    const unsigned long long p=static_cast<unsigned long long>(optixGetPayload_0())
        | (static_cast<unsigned long long>(optixGetPayload_1())<<32);
    return reinterpret_cast<QueryState*>(p);
}
}

extern "C" __global__ void __raygen__patch_query()
{
    const unsigned index=optixGetLaunchIndex().x;
    if(index>=patch_query_params.direction_count) return;
    QueryState state{};
    state.collector.writing=patch_query_params.write_pass!=0;
    if(state.collector.writing)
    {
        const auto begin=patch_query_params.offsets[index],end=patch_query_params.offsets[index+1];
        if(end<begin || end-begin>0xffffffffull)
        {
            rainbow::PatchQuerySummary failure{};
            failure.flags=rainbow::PatchQueryInvalidData;
            patch_query_params.summaries[index]=failure;
            return;
        }
        state.collector.capacity=static_cast<std::uint32_t>(end-begin);
        if(end!=begin) state.collector.storage=patch_query_params.hits+begin;
    }
    const auto direction=patch_query_params.directions[index];
    const float norm2=direction.dot(direction);
    if(!direction.is_finite() || !(norm2>0.99f && norm2<1.01f))
        state.collector.flags|=rainbow::PatchQueryInvalidData;
    else if(patch_query_params.traversable!=0)
    {
        const auto address=reinterpret_cast<unsigned long long>(&state);
        unsigned low=static_cast<unsigned>(address),high=static_cast<unsigned>(address>>32);
        // All-hit: AH を有効にし，closest-hit は使わない．front/back culling もしない．
        // corners は単位方向空間なので，正の交点は t<4 に入る．
        optixTrace(static_cast<OptixTraversableHandle>(patch_query_params.traversable),
            make_float3(0,0,0),make_float3(direction.x,direction.y,direction.z),
            0.0f,4.0f,0.0f,255,OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT,
            0,1,0,low,high);
    }
    if(state.collector.writing)
        patch_query_params.summaries[index]=state.collector.finish(
            patch_query_params.patches,patch_query_params.patch_count);
    else
    {
        rainbow::PatchQuerySummary count{};
        count.candidate_count=state.collector.count;
        count.flags=state.collector.flags;
        count.first_problem_patch_id=state.collector.first_problem_patch_id;
        patch_query_params.summaries[index]=count;
    }
}

extern "C" __global__ void __intersection__patch_query()
{
    auto* state=query_state();
    const unsigned primitive=optixGetPrimitiveIndex();
    if(primitive>=patch_query_params.patch_count)
    {state->collector.flags|=rainbow::PatchQueryInvalidData;return;}
    const auto& patch=patch_query_params.patches[primitive];
    rainbow::BilinearPatchGeometry geometry{};
    for(unsigned i=0;i<4;++i)
    {
        if(patch.vertex_indices[i]>=patch_query_params.vertex_count)
        {state->collector.flags|=rainbow::PatchQueryInvalidData;return;}
        geometry.corners[i]=patch_query_params.vertices[patch.vertex_indices[i]].direction_drop;
    }
    const auto w=optixGetObjectRayDirection();
    state->pending=rainbow::BilinearPatchIntersector::intersect(geometry,{w.x,w.y,w.z});
    // 解けない候補を，silent miss として消さない．これは OR のみなので訪問順序に依らない．
    state->collector.note(state->pending.flags,patch.patch_id);
    if(state->pending.count==0) return;
    state->compact_index=primitive;
    state->consumed=false;
    float report_t=state->pending.hits[0].t;
    if(state->pending.count==2 && state->pending.hits[1].t<report_t)
        report_t=state->pending.hits[1].t;
    // 1 primitive につき1回 report し，AH で最大2個の根を両方取り込む．
    // REQUIRE_SINGLE_ANYHIT_CALL と組み合わせ，BVH 内の空間分割による二重収集を避ける．
    static_cast<void>(optixReportIntersection(report_t,0));
    // AH は ignore を呼ぶので report の戻り値 false は正常．consumed で呼出しを確認する．
    if(!state->consumed) state->collector.flags|=rainbow::PatchQueryInvalidData;
}

extern "C" __global__ void __anyhit__patch_query()
{
    auto* state=query_state();
    const auto primitive=state->compact_index;
    state->collector.consume(patch_query_params.patches[primitive],primitive,state->pending);
    state->consumed=true;
    // hit を記録すると tmax が縮まり遠い patch を失う．全件取得のため必ず ignore．
    optixIgnoreIntersection();
}

extern "C" __global__ void __miss__patch_query()
{
    // 全候補を ignore するため，ヒットが何件あっても最終的に miss に来る．
    // collector をゼロに戻してはいけない．
}
