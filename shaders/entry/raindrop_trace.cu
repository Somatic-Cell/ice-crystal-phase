#include <optix.h>
#include <rainbow/trace_launch_params.hpp>

extern "C"
{
    __constant__ rainbow::TraceLaunchParams params;
}

extern "C" __global__
void __raygen__raindrop_trace()
{
    const unsigned int index = optixGetLaunchIndex().x;
    if(index < params.count)
    {
        params.output[index] =
            index ^ 0xa5a5a5a5u;
    }
}

extern "C" __global__
void __intersection__raindrop()
{

}

extern "C" __global__
void __closesthit__raindrop()
{

}

extern "C" __global__
void __miss__raindrop()
{

}