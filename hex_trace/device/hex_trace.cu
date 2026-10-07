#include <ice_crystal/hex_trace_core.hpp>
using namespace iceCrystal;

// ABI fingerprint includes structure sizes, so a stale module is rejected even
// before the first expensive tracing pass.
extern "C" __global__ void ice_hex_trace_abi_v1(std::uint64_t* words)
{
    if(blockIdx.x||threadIdx.x)return;
    words[0]=UINT64_C(0x4943454845580001);
    words[1]=sizeof(TraceSettings);words[2]=sizeof(RayAudit);words[3]=sizeof(OutgoingSample);
}
extern "C" __global__ void ice_hex_trace_count(TraceSettings s,std::uint64_t first,
    std::uint64_t count,RayAudit* audits)
{
    for(std::uint64_t i=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;i<count;
        i+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        CountSink sink{};audits[i]=trace_sample(s,first+i,sink);
    }
}
// Ordered integer prefix, exact unless overflow is reported. The array is only
// one entry per INCIDENT ray, not per outgoing direction or phase-function cell.
// Kept simple for M1; no change of estimator is involved.
extern "C" __global__ void ice_hex_trace_offsets(const RayAudit* audits,std::uint64_t n,
    std::uint64_t* offsets,std::uint64_t* total)
{
    if(blockIdx.x||threadIdx.x)return;
    offsets[0]=0;std::uint64_t sum=0;
    for(std::uint64_t i=0;i<n;++i)
    {
        const auto c=audits[i].output_count;
        if(c>~std::uint64_t{0}-sum){*total=~std::uint64_t{0};return;}
        sum+=c;offsets[i+1]=sum;
    }
    *total=sum;
}
extern "C" __global__ void ice_hex_trace_write(TraceSettings s,std::uint64_t first,
    std::uint64_t count,const std::uint64_t* offsets,OutgoingSample* output,RayAudit* audits)
{
    for(std::uint64_t i=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;i<count;
        i+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        const auto n=offsets[i+1]-offsets[i];
        BufferSink sink{n?output+offsets[i]:nullptr,n,0};
        const auto replay=trace_sample(s,first+i,sink);
        if(!replay_equal(audits[i],replay)||sink.size!=n)
            audits[i].status=TraceStatus::ReplayMismatch;
    }
}
