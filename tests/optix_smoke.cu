// 旧 raindrop_trace.cu の M2 raygen．既存テストの ABI/期待値/entry 名を維持する．
#include <optix.h>
#include <rainbow/trace_launch_params.hpp>
extern "C" { __constant__ rainbow::TraceLaunchParams params; }
extern "C" __global__ void __raygen__raindrop_trace()
{
    const unsigned index=optixGetLaunchIndex().x;
    if(index<params.count) params.output[index]=index^0xa5a5a5a5u;
}
