#include <rainbow/phase_cdf_data.hpp>
#include <cmath>

namespace
{
using namespace rainbow;
constexpr unsigned threads = 256;

__device__ PhaseCdfReduction empty_reduction()
{
    PhaseCdfReduction s{};
    s.first_problem = ~std::uint64_t{0};
    return s;
}
__device__ PhaseCdfReduction combine(PhaseCdfReduction a, const PhaseCdfReduction& b)
{
    a.maximum = fmax(a.maximum, b.maximum);
    a.l1_error += b.l1_error;
    a.lost_mass += b.lost_mass;
    a.maximum_cell_error = fmax(a.maximum_cell_error, b.maximum_cell_error);
    a.invalid_values += b.invalid_values; a.incomplete += b.incomplete;
    a.underresolved += b.underresolved; a.malformed += b.malformed;
    a.lost_cells += b.lost_cells; a.weight_underflow += b.weight_underflow;
    a.first_problem = a.first_problem < b.first_problem ? a.first_problem : b.first_problem;
    return a;
}
__device__ void block_reduce(PhaseCdfReduction value, PhaseCdfReduction* destination)
{
    __shared__ PhaseCdfReduction shared[threads];
    shared[threadIdx.x] = value;
    __syncthreads();
    for(unsigned step=threads/2; step; step/=2)
    {
        if(threadIdx.x < step)
            shared[threadIdx.x] = combine(shared[threadIdx.x], shared[threadIdx.x+step]);
        __syncthreads();
    }
    if(threadIdx.x == 0) *destination = shared[0];
}
__device__ double density(PhaseDensityView v, std::uint64_t i,
                          bool& incomplete, bool& underresolved)
{
    incomplete = false; underresolved = false;
    if(v.stage == PhaseDensityStage::Scalar) return v.scalar[i];
    incomplete = !v.optical[i].known_hits_complete();
    if(v.stage == PhaseDensityStage::Incoherent)
        return v.optical[i].regular_partial_incoherent_s + v.optical[i].regular_partial_incoherent_p;
    if(v.stage == PhaseDensityStage::Path)
        return v.optical[i].regular_partial_path_s + v.optical[i].regular_partial_path_p;
    incomplete = incomplete || !v.focal[i].valid();
    if(v.stage == PhaseDensityStage::Focal)
        return v.focal[i].intensity_s + v.focal[i].intensity_p;
    incomplete = incomplete || !v.diffraction[i].valid();
    underresolved = (v.diffraction[i].flags & DiffractionUnderresolved) != 0;
    return v.diffraction[i].intensity_s + v.diffraction[i].intensity_p;
}
__device__ bool finite_field(const Field32& f)
{
    return isfinite(f.x.real) && isfinite(f.x.imag) && isfinite(f.y.real) && isfinite(f.y.imag);
}


}

extern "C" __global__ void rainbow_phase_cdf_abi_v1() {}

extern "C" __global__ void phase_validate_trace(rainbow::TraceValidationParams p)
{
    PhaseCdfReduction local=empty_reduction();
    for(std::uint64_t i=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        i<p.count; i+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        const auto& v=p.vertices[i];
        bool bad=v.status!=VertexStatus::Valid && v.status!=VertexStatus::Miss &&
                 v.status!=VertexStatus::TotalInternalReflection;
        if(v.status==VertexStatus::Valid)
            bad=bad || !v.position_drop.is_finite() || !v.direction_drop.is_finite() ||
                !v.basis_x.is_finite() || !v.optical_cycles.is_valid() ||
                !finite_field(v.field) || !finite_field(p.second[i]);
        if(bad) { ++local.incomplete; local.first_problem=(local.first_problem<i?local.first_problem:i); }
    }
    block_reduce(local, p.partials+blockIdx.x);
}

extern "C" __global__ void phase_reduce_summary(const rainbow::PhaseCdfReduction* input,
    std::uint32_t count, rainbow::PhaseCdfReduction* output)
{
    PhaseCdfReduction local=empty_reduction();
    for(std::uint32_t i=threadIdx.x; i<count; i+=blockDim.x) local=combine(local,input[i]);
    block_reduce(local,output);
}

extern "C" __global__ void phase_check_density(rainbow::PhaseCdfBuildParams p)
{
    PhaseCdfReduction local=empty_reduction();
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<p.input.count; k+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        bool incomplete=false, under=false;
        const double d=density(p.input,k,incomplete,under);
        const bool bad=!isfinite(d) || d<0;
        local.invalid_values+=bad; local.incomplete+=incomplete; local.underresolved+=under;
        if(bad || incomplete || under)
            local.first_problem=local.first_problem<k?local.first_problem:k;
        if(!bad) local.maximum=fmax(local.maximum,d);
    }
    block_reduce(local,p.partials+blockIdx.x);
}

// Failure-path inspection only. The density and completeness predicates are
// exactly those of phase_check_density. Only the meaning of first_problem is
// different: an underresolution warning alone must not hide a hard failure.
// No optical value, validity flag, CDF element or acceptance policy is modified.
extern "C" __global__ void phase_check_density_failures_only(rainbow::PhaseCdfBuildParams p)
{
    PhaseCdfReduction local=empty_reduction();
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<p.input.count; k+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        bool incomplete=false, under=false;
        const double d=density(p.input,k,incomplete,under);
        const bool bad=!isfinite(d) || d<0;
        local.invalid_values+=bad;
        local.incomplete+=incomplete;
        local.underresolved+=under;
        if(bad || incomplete)
            local.first_problem=local.first_problem<k?local.first_problem:k;
        if(!bad) local.maximum=fmax(local.maximum,d);
    }
    block_reduce(local,p.partials+blockIdx.x);
}

// Transpose row-major intensity into phi-major conditional storage.
// dOmega = (4*pi/Nphi)*du; the common constant cancels during normalization.
extern "C" __global__ void phase_make_weights(rainbow::PhaseCdfBuildParams p)
{
    __shared__ double tile[32][33];
    const unsigned col=blockIdx.x*32+threadIdx.x;
    const unsigned row=blockIdx.y*32+threadIdx.y;
    for(unsigned dy=0; dy<32; dy+=8)
    {
        double w=0;
        if(col<p.phi_count && row+dy<p.theta_count)
        {
            bool a=false,b=false;
            const double d=density(p.input,std::uint64_t(row+dy)*p.phi_count+col,a,b);
            w=(d/p.input_summary->maximum)*(p.u_edges[row+dy+1]-p.u_edges[row+dy]);
        }
        tile[threadIdx.y+dy][threadIdx.x]=w;
    }
    __syncthreads();
    const unsigned out_phi=blockIdx.x*32+threadIdx.y;
    const unsigned out_theta=blockIdx.y*32+threadIdx.x;
    for(unsigned dy=0; dy<32; dy+=8)
        if(out_phi+dy<p.phi_count && out_theta<p.theta_count)
            p.theta_cdf[std::uint64_t(out_phi+dy)*(p.theta_count+1ull)+out_theta+1ull]=
                tile[threadIdx.x][threadIdx.y+dy];
}
// Columns are independent and evaluated in parallel. Within each column use
// ordered nonnegative binary64 additions: exact zero weights give exact CDF
// plateaus and rounding cannot make a prefix decrease. A parallel re-associated
// scan can violate both properties in extreme high-dynamic-range data. Error
// relative to the original masses is measured by phase_audit_cdf; never repaired.
extern "C" __global__ void phase_scan_conditionals(rainbow::PhaseCdfBuildParams p)
{
    const std::uint64_t j=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(j>=p.phi_count)return;
    const std::uint64_t base=j*(p.theta_count+1ull);
    p.theta_cdf[base]=0;
    double sum=0;
    for(std::uint32_t i=0;i<p.theta_count;++i)
    {
        sum+=p.theta_cdf[base+i+1];
        p.theta_cdf[base+i+1]=sum;
    }
    p.column_sums[j]=sum;
}
// Marginal contains only Nphi numbers (not Ntheta*Nphi). A single ordered scan
// retains the same plateau/monotonicity property, without host synchronization.
extern "C" __global__ void phase_scan_marginal(rainbow::PhaseCdfBuildParams p)
{
    if(blockIdx.x!=0 || threadIdx.x!=0)return;
    p.phi_cdf[0]=0;
    double sum=0;
    for(std::uint32_t j=0;j<p.phi_count;++j)
    {
        sum+=p.column_sums[j];
        p.phi_cdf[j+1]=sum;
    }
    *p.total=sum;
}
extern "C" __global__ void phase_normalize(rainbow::PhaseCdfBuildParams p)
{
    const std::uint64_t total=std::uint64_t(p.phi_count)*(p.theta_count+1ull);
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<total; k+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        const auto j=k/(p.theta_count+1ull),i=k%(p.theta_count+1ull);
        if(i==0) p.theta_cdf[k]=0;
        else if(i==p.theta_count) p.theta_cdf[k]=1;
        else p.theta_cdf[k]=p.column_sums[j]>0?p.theta_cdf[k]/p.column_sums[j]:p.u_edges[i];
    }
    for(std::uint64_t j=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        j<=p.phi_count; j+=std::uint64_t(gridDim.x)*blockDim.x)
        p.phi_cdf[j]=j==0?0:(j==p.phi_count?1:p.phi_cdf[j]/(*p.total));
}
extern "C" __global__ void phase_audit_cdf(rainbow::PhaseCdfBuildParams p)
{
    PhaseCdfReduction local=empty_reduction();
    const double total=*p.total;
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<p.input.count; k+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        const auto i=k/p.phi_count,j=k%p.phi_count;
        const auto base=j*(p.theta_count+1ull)+i;
        const double a0=p.phi_cdf[j],a1=p.phi_cdf[j+1];
        const double b0=p.theta_cdf[base],b1=p.theta_cdf[base+1];
        const bool bad=!(isfinite(a0)&&isfinite(a1)&&isfinite(b0)&&isfinite(b1)&&
            a0>=0&&a1<=1&&a1>=a0&&b0>=0&&b1<=1&&b1>=b0&&total>0&&isfinite(total)) ||
            (j==0&&a0!=0) || (j+1==p.phi_count&&a1!=1) ||
            (i==0&&b0!=0) || (i+1==p.theta_count&&b1!=1);
        local.malformed+=bad;
        if(bad) { local.first_problem=local.first_problem<k?local.first_problem:k; continue; }
        bool unused_a=false,unused_b=false;
        const double d=density(p.input,k,unused_a,unused_b);
        const double w=(d/p.input_summary->maximum)*(p.u_edges[i+1]-p.u_edges[i]);
        const double original=w/total, reconstructed=(a1-a0)*(b1-b0);
        const bool under=d>0 && (!(w>0)||!(original>0));
        local.weight_underflow+=under;
        if(under) local.first_problem=local.first_problem<k?local.first_problem:k;
        const double error=fabs(reconstructed-original);
        local.l1_error+=error; local.maximum_cell_error=fmax(local.maximum_cell_error,error);
        if(original>0 && reconstructed==0) { ++local.lost_cells; local.lost_mass+=original; }
    }
    block_reduce(local,p.partials+blockIdx.x);
}

// Storage v2: spherical prefilter, conservative aggregation, HG first moment.
#include <rainbow/phase_storage_math.hpp>
extern "C" __global__ void rainbow_phase_storage_abi_v2() {}
namespace
{
__device__ void reduce_moments(rainbow::PhaseMomentSum a, rainbow::PhaseMomentSum* out)
{
    __shared__ rainbow::PhaseMomentSum shared[256];
    shared[threadIdx.x]=a;
    __syncthreads();
    for(unsigned s=128;s;s/=2)
    {
        if(threadIdx.x<s)
        {
            shared[threadIdx.x].mass+=shared[threadIdx.x+s].mass;
            shared[threadIdx.x].axial+=shared[threadIdx.x+s].axial;
            shared[threadIdx.x].tv+=shared[threadIdx.x+s].tv;
            shared[threadIdx.x].weight_underflow+=shared[threadIdx.x+s].weight_underflow;
        }
        __syncthreads();
    }
    if(threadIdx.x==0)*out=shared[0];
}
}
extern "C" __global__ void phase_storage_reduce_moments(
    const rainbow::PhaseMomentSum* input, std::uint32_t count, rainbow::PhaseMomentSum* output)
{
    rainbow::phase_storage_math::Sum mass{},axial{},tv{};
    std::uint64_t underflow=0;
    for(std::uint32_t k=threadIdx.x;k<count;k+=blockDim.x)
    {mass.add(input[k].mass);axial.add(input[k].axial);tv.add(input[k].tv);underflow+=input[k].weight_underflow;}
    reduce_moments({mass.value(),axial.value(),tv.value(),underflow},output);
}
extern "C" __global__ void phase_storage_extract(rainbow::PhaseStorageParams p)
{
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<p.input.count;k+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        bool a=false,b=false;
        p.output[k]=density(p.input,k,a,b)/p.divisor;
    }
}
extern "C" __global__ void phase_storage_gaussian(rainbow::PhaseStorageParams p)
{
    const rainbow::phase_storage_math::Grid g{p.fine_values,p.fine_edges,p.nt,p.np};
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<std::uint64_t(p.nt)*p.np;k+=std::uint64_t(gridDim.x)*blockDim.x)
        p.output[k]=rainbow::phase_storage_math::gaussian_at(
            g,std::uint32_t(k/p.np),std::uint32_t(k%p.np),p.sigma_rad,p.support_sigma);
}
extern "C" __global__ void phase_storage_coarsen(rainbow::PhaseStorageParams p)
{
    const rainbow::phase_storage_math::Grid g{p.fine_values,p.fine_edges,p.nt,p.np};
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<std::uint64_t(p.out_nt)*p.out_np;k+=std::uint64_t(gridDim.x)*blockDim.x)
        p.output[k]=rainbow::phase_storage_math::coarsen_at(
            g,std::uint32_t(k/p.out_np),std::uint32_t(k%p.out_np),p.out_nt,p.out_np);
}
extern "C" __global__ void phase_storage_measure(rainbow::PhaseStorageParams p)
{
    rainbow::phase_storage_math::Sum mass{},axial{};
    std::uint64_t underflow=0;
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<p.input.count;k+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        bool a=false,b=false;
        const auto i=std::uint32_t(k/p.np);
        const double d=density(p.input,k,a,b);
        const double weight=(d/p.divisor)*(p.fine_edges[i+1]-p.fine_edges[i]);
        underflow+=d>0 && !(weight>0);
        const double w=weight/double(p.np);
        mass.add(w);axial.add(w*rainbow::phase_storage_math::mean_cosine(p.fine_edges[i],p.fine_edges[i+1]));
    }
    reduce_moments({mass.value(),axial.value(),0,underflow},p.moment_partials+blockIdx.x);
}
extern "C" __global__ void phase_storage_measure_tv(rainbow::PhaseStorageParams p)
{
    rainbow::phase_storage_math::Sum tv{};
    const auto rt=p.nt/p.out_nt,rp=p.np/p.out_np;
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<std::uint64_t(p.nt)*p.np;k+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        const auto i=std::uint32_t(k/p.np),j=std::uint32_t(k%p.np);
        const double a=p.fine_values[k]/p.fine_mass;
        const double b=p.coarse_values[std::uint64_t(i/rt)*p.out_np+j/rp]/p.coarse_mass;
        tv.add(0.5*::fabs(a-b)*(p.fine_edges[i+1]-p.fine_edges[i])/double(p.np));
    }
    reduce_moments({0,0,tv.value(),0},p.moment_partials+blockIdx.x);
}
extern "C" __global__ void phase_storage_measure_cdf(rainbow::PhaseStorageParams p)
{
    rainbow::phase_storage_math::Sum mass{},axial{};
    for(std::uint64_t k=std::uint64_t(blockIdx.x)*blockDim.x+threadIdx.x;
        k<std::uint64_t(p.out_nt)*p.out_np;k+=std::uint64_t(gridDim.x)*blockDim.x)
    {
        const auto i=std::uint32_t(k/p.out_np),j=std::uint32_t(k%p.out_np);
        const auto b=std::uint64_t(j)*(p.out_nt+1ull)+i;
        const double w=(p.phi_cdf[j+1]-p.phi_cdf[j])*(p.theta_cdf[b+1]-p.theta_cdf[b]);
        mass.add(w);axial.add(w*rainbow::phase_storage_math::mean_cosine(p.coarse_edges[i],p.coarse_edges[i+1]));
    }
    reduce_moments({mass.value(),axial.value(),0,0},p.moment_partials+blockIdx.x);
}
