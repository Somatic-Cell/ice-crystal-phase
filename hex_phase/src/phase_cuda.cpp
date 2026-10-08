#include <ice_crystal/phase_cuda.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/cuda_host_upload.hpp>
#include <algorithm>
#include <stdexcept>

namespace iceCrystal
{
PhaseMassAccumulatorCuda::PhaseMassAccumulatorCuda(const rainbow::CudaContext& c,PhaseGrid g)
    :grid_(std::move(g)),context_(c),module_(c),edges_(c),mass_(c),density_(c),statistics_(c),partials_(c)
{
    validate_phase_grid(grid_);const auto n=std::size_t(grid_.nt)*grid_.np;
    edges_.allocate(grid_.edges.size());mass_.allocate(n);density_.allocate(n);statistics_.allocate(1);
    rainbow::upload_host_to_device_sync(context_, edges_.address(), grid_.edges.data(), edges_.byte_size());
    mass_.zero_byte_async(context_.stream());statistics_.zero_byte_async(context_.stream());pending_=true;sync();
}
PhaseMassAccumulatorCuda::~PhaseMassAccumulatorCuda() noexcept
{
    if(pending_)
    {
        static_cast<void>(rainbow::detail::report_cuda_cleanup_result(cuCtxSetCurrent(context_.handle()),"phase accumulator cleanup context"));
        static_cast<void>(rainbow::detail::report_cuda_cleanup_result(cuStreamSynchronize(context_.stream()),"phase accumulator cleanup stream"));
    }
}
void PhaseMassAccumulatorCuda::sync()
{context_.make_current();RAINBOW_CUDA_CHECK(cuStreamSynchronize(context_.stream()));pending_=false;}
void PhaseMassAccumulatorCuda::launch(const char* name,void** args,unsigned blocks,unsigned threads)
{
    context_.make_current();const auto f=module_.find_function(name);pending_=true;
    RAINBOW_CUDA_CHECK(cuLaunchKernel(f,blocks,1,1,threads,1,1,0,context_.stream(),args,nullptr));
}
void PhaseMassAccumulatorCuda::load_module(const std::filesystem::path& path)
{
    module_.load_fatbin(path);static_cast<void>(module_.find_function("ice_phase_accumulate_abi_m2_v1"));
    static_cast<void>(module_.find_function("ice_phase_accumulate"));static_cast<void>(module_.find_function("ice_phase_finish"));ready_=true;
}
void PhaseMassAccumulatorCuda::add(const rainbow::DeviceBuffer<OutgoingSample>& samples,const Rotation& r,double scale)
{
    if(!ready_||finished_)throw std::logic_error("Phase accumulator not ready / already finished.");
    if(!r.valid()||!(scale>0)||!finite_value(scale))throw std::invalid_argument("Invalid orientation / area weight.");
    if(samples.is_empty())return;
    PhaseAccumulateParams p{};p.grid=grid_.view();p.grid.u_edges=edges_.data();p.rotation=r;
    p.outgoing=reinterpret_cast<const OutgoingSample*>(samples.address());p.count=samples.element_count();
    p.cell_mass=mass_.data();p.statistics=statistics_.data();p.area_weight_per_ray=scale;
    void* args[]{&p};launch("ice_phase_accumulate",args,static_cast<unsigned>(std::min<std::uint64_t>(4096,(p.count+255)/256)));
    // Synchronous borrowing: the M1 tracer may release its previous output on the next batch.
    sync();
}
void PhaseMassAccumulatorCuda::finish()
{
    if(!ready_||finished_)throw std::logic_error("Phase accumulator not ready / already finished.");
    sync();statistics_.download({&host_statistics_,1});
    if(host_statistics_.invalid_samples||!finite_value(host_statistics_.point_mass)||!(host_statistics_.point_mass>0)||
       !finite_value(host_statistics_.point_axial)||!finite_value(host_statistics_.forward_mass)||
       !finite_value(host_statistics_.backward_mass)||!finite_value(host_statistics_.absorbed_addend_mass))
        throw std::runtime_error("GPU histogram rejected invalid samples / nonfinite accumulated powers.");
    const auto n=std::uint64_t(grid_.nt)*grid_.np;const auto blocks=static_cast<unsigned>(std::min<std::uint64_t>(4096,(n+255)/256));
    partials_.allocate(blocks);PhaseFinishParams p{};p.grid=grid_.view();p.grid.u_edges=edges_.data();
    p.cell_mass=mass_.data();p.statistics=statistics_.data();p.density=density_.data();p.partials=partials_.data();
    void* args[]{&p};launch("ice_phase_finish",args,blocks);sync();
    std::vector<HistogramMoments> values(blocks);partials_.download(values);Sum mass{},axial{};
    for(const auto& a:values){mass.add(a.mass);axial.add(a.axial);}
    host_moments_={mass.value,axial.value};
    if(!(mass.value>0)||!finite_value(mass.value)||!finite_value(axial.value))throw std::runtime_error("GPU histogram integral invalid.");
    finished_=true;
}
HistogramStats PhaseMassAccumulatorCuda::statistics() const
{if(!finished_)throw std::logic_error("Histogram not finalized.");return host_statistics_;}
HistogramMoments PhaseMassAccumulatorCuda::moments() const
{if(!finished_)throw std::logic_error("Histogram not finalized.");return host_moments_;}
const rainbow::DeviceBuffer<double>& PhaseMassAccumulatorCuda::density() const
{if(!finished_)throw std::logic_error("Histogram not finalized.");return density_;}
std::vector<double> PhaseMassAccumulatorCuda::download_density() const
{const auto& d=density();std::vector<double> v(d.element_count());d.download(v);return v;}
} // namespace iceCrystal
