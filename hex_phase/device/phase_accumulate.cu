#include <ice_crystal/phase_math.hpp>

namespace
{
using namespace iceCrystal;
__device__ void increment(std::uint64_t* p)
{atomicAdd(reinterpret_cast<unsigned long long*>(p),1ull);}
}
extern "C" __global__ void ice_phase_accumulate_abi_m2_v1() {}
extern "C" __global__ void ice_phase_accumulate(iceCrystal::PhaseAccumulateParams p)
{
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<p.count;k+=std::uint64_t(blockDim.x)*gridDim.x)
    {
        const auto& s=p.outgoing[k];const double w=p.area_weight_per_ray*s.power_fraction;
        const auto b=phase_bin(p.grid,p.rotation.apply(s.direction));
        if(!b.valid||!(s.power_fraction>0)||!finite_value(s.power_fraction)||!(w>0)||!finite_value(w))
        {increment(&p.statistics->invalid_samples);continue;}
        double* target=b.pole>0?&p.statistics->forward_mass:
            (b.pole<0?&p.statistics->backward_mass:p.cell_mass+b.index);
        const double old=atomicAdd(target,w);
        if(old+w==old)
        {increment(&p.statistics->absorbed_addends);atomicAdd(&p.statistics->absorbed_addend_mass,w);}
        atomicAdd(&p.statistics->point_mass,w);atomicAdd(&p.statistics->point_axial,w*b.mu);
        increment(&p.statistics->samples);
        if(b.boundary_snapped)increment(&p.statistics->boundary_snapped_samples);
    }
}
extern "C" __global__ void ice_phase_finish(iceCrystal::PhaseFinishParams p)
{
    __shared__ double shared_mass[256],shared_axial[256];
    Sum mass{},axial{};const auto n=std::uint64_t(p.grid.nt)*p.grid.np;
    const double forward=p.statistics->forward_mass/p.grid.np,backward=p.statistics->backward_mass/p.grid.np;
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;k<n;k+=std::uint64_t(blockDim.x)*gridDim.x)
    {
        const auto i=k/p.grid.np;double w=p.cell_mass[k];
        if(i==0)w+=forward;
        if(i+1==p.grid.nt)w+=backward;
        const double d=w/((4*pi/p.grid.np)*(p.grid.u_edges[i+1]-p.grid.u_edges[i]));
        // Preserve a failure as NaN for the existing PhaseCdf validity audit.
        p.density[k]=(!(d>=0)||!finite_value(d)||(w>0&&!(d>0))||
            (p.statistics->forward_mass>0&&!(forward>0))||
            (p.statistics->backward_mass>0&&!(backward>0)))?::nan(""):d;
        mass.add(w);axial.add(w*(1-p.grid.u_edges[i]-p.grid.u_edges[i+1]));
    }
    shared_mass[threadIdx.x]=mass.value;shared_axial[threadIdx.x]=axial.value;__syncthreads();
    for(unsigned step=128;step;step/=2)
    {
        if(threadIdx.x<step){shared_mass[threadIdx.x]+=shared_mass[threadIdx.x+step];shared_axial[threadIdx.x]+=shared_axial[threadIdx.x+step];}
        __syncthreads();
    }
    if(threadIdx.x==0)p.partials[blockIdx.x]={shared_mass[0],shared_axial[0]};
}
