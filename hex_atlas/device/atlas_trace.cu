#include <ice_crystal/atlas.hpp>
extern "C" __global__ void ice_atlas_abi_m3_v1() {}
extern "C" __global__ void atlas_trace_area(iceCrystal::atlas::BatchParams p)
{
    const auto i=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i>=p.count)return;
    const auto q=p.rays[i];
    iceCrystal::atlas::MassSink sink{p.grid,p.mass,p.count,i,q.projected_weight_mm2};
    p.audits[i]=iceCrystal::trace_from_entry(p.trace,p.first+i,q.entry,sink);
}
extern "C" __global__ void atlas_reduce_ordered(const iceCrystal::Sum* partial,
    iceCrystal::Sum* total,std::uint64_t cells,std::uint64_t count)
{
    const auto bin=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(bin>=cells)return;
    auto sum=total[bin];
    for(std::uint64_t i=0;i<count;++i)
    {
        const double w=partial[bin*count+i].value;
        if(w>0)sum.add(w);
    }
    total[bin]=sum;
}
