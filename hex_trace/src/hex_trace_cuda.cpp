#include <ice_crystal/hex_trace_cuda.hpp>
#include <rainbow/cuda_error.hpp>
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>

namespace iceCrystal
{
HexPrismTracer::HexPrismTracer(const rainbow::CudaContext& c)
    :context_(c),module_(c),audits_(c),offsets_(c),total_(c),outgoing_(c) {}
HexPrismTracer::~HexPrismTracer() noexcept
{
    if(pending_)
    {
        static_cast<void>(rainbow::detail::report_cuda_cleanup_result(
            cuCtxSetCurrent(context_.handle()),"hex tracer cleanup context"));
        static_cast<void>(rainbow::detail::report_cuda_cleanup_result(
            cuStreamSynchronize(context_.stream()),"hex tracer cleanup stream"));
    }
}
void HexPrismTracer::synchronize()
{
    context_.make_current();RAINBOW_CUDA_CHECK(cuStreamSynchronize(context_.stream()));pending_=false;
}
void HexPrismTracer::launch(const char* name,void** args,unsigned blocks)
{
    context_.make_current();const auto f=module_.find_function(name);pending_=true;
    RAINBOW_CUDA_CHECK(cuLaunchKernel(f,blocks,1,1,128,1,1,0,context_.stream(),args,nullptr));
}
void HexPrismTracer::load_module(const std::filesystem::path& path)
{
    module_ready_=false;
    module_.load_fatbin(path);
    static_cast<void>(module_.find_function("ice_hex_trace_count"));
    static_cast<void>(module_.find_function("ice_hex_trace_offsets"));
    static_cast<void>(module_.find_function("ice_hex_trace_write"));
    // Use member-owned memory so even a failed launch/sync has a safe lifetime.
    total_.allocate(4);auto* ptr=total_.data();void* args[]={&ptr};
    launch("ice_hex_trace_abi_v1",args,1);synchronize();
    std::array<std::uint64_t,4> abi{};total_.download(abi);total_.close();
    const std::array<std::uint64_t,4> expected{
        UINT64_C(0x4943454845580001),sizeof(TraceSettings),sizeof(RayAudit),sizeof(OutgoingSample)};
    if(abi!=expected)throw std::runtime_error("Hex trace fatbin ABI mismatch; rebuild/stage the module.");
    module_ready_=true;
}
void HexPrismTracer::trace(const TraceSettings& config,std::uint64_t first,std::size_t count,std::size_t max_records)
{
    has_result_=false;validate_batch_range(config,first,count);
    if(!module_ready_)throw std::logic_error("Hex trace module not loaded.");
    synchronize();audits_.close();offsets_.close();total_.close();outgoing_.close();
    audits_.allocate(count);offsets_.allocate(count+1);total_.allocate(1);
    auto s=config;const auto n=static_cast<std::uint64_t>(count);
    auto* a=audits_.data();auto* o=offsets_.data();auto* t=total_.data();
    auto first_arg=first,n_arg=n;
    void* count_args[]={&s,&first_arg,&n_arg,&a};
    const auto blocks=static_cast<unsigned>((std::min)(std::uint64_t{4096},n/128+(n%128!=0)));
    launch("ice_hex_trace_count",count_args,blocks);
    void* prefix_args[]={&a,&n_arg,&o,&t};launch("ice_hex_trace_offsets",prefix_args,1);synchronize();
    std::uint64_t output_count=0;total_.download({&output_count,1});
    if(output_count>max_records || output_count>(std::numeric_limits<std::size_t>::max)()/sizeof(OutgoingSample))
        throw std::length_error("Hex output memory budget/size exceeded; no automatic sample reduction.");
    outgoing_.allocate(static_cast<std::size_t>(output_count));
    pending_=true;outgoing_.zero_byte_async(context_.stream());
    auto* output=outgoing_.data();
    void* write_args[]={&s,&first_arg,&n_arg,&o,&output,&a};
    launch("ice_hex_trace_write",write_args,blocks);synchronize();
    count_=count;first_=first;has_result_=true;
}
std::vector<RayAudit> HexPrismTracer::download_audits() const
{
    if(!has_result_)throw std::logic_error("No completed hex trace buffers.");
    std::vector<RayAudit> a(count_);audits_.download(a);return a;
}
TraceBatch HexPrismTracer::download() const
{
    if(!has_result_)throw std::logic_error("No completed hex trace buffers.");
    TraceBatch b{};b.first_sample=first_;b.audits=download_audits();
    b.offsets.resize(count_+1);offsets_.download(b.offsets);
    b.outgoing.resize(outgoing_.element_count());outgoing_.download(b.outgoing);return b;
}
} // namespace iceCrystal
