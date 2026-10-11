#include <ice_crystal/atlas_cuda.hpp>
#include <rainbow/cuda_host_upload.hpp>
#include <algorithm>
#include <limits>
#include <stdexcept>
namespace iceCrystal::atlas
{
CudaTracer::CudaTracer(const rainbow::CudaContext& c,const std::filesystem::path& path)
    :context_(c),module_(c),rays_(c),audits_(c),partial_(c),mass_(c),edges_(c)
{
    module_.load_fatbin(path);static_cast<void>(module_.find_function("ice_atlas_abi_m3_v1"));
}
CudaTracer::~CudaTracer() noexcept
{
    if(pending_)
    {
        static_cast<void>(rainbow::detail::report_cuda_cleanup_result(cuCtxSetCurrent(context_.handle()),"atlas cleanup context"));
        static_cast<void>(rainbow::detail::report_cuda_cleanup_result(cuStreamSynchronize(context_.stream()),"atlas cleanup sync"));
    }
}
void CudaTracer::sync()
{context_.make_current();RAINBOW_CUDA_CHECK(cuStreamSynchronize(context_.stream()));pending_=false;}
void CudaTracer::launch(const char* name,void** args,unsigned blocks)
{
    context_.make_current();pending_=true;
    RAINBOW_CUDA_CHECK(cuLaunchKernel(module_.find_function(name),blocks,1,1,128,1,1,0,context_.stream(),args,nullptr));
}
Result CudaTracer::trace(const Settings& s,const Node& n)
{
    validate_node(s,n);sync();rays_.close();audits_.close();partial_.close();mass_.close();edges_.close();
    const auto batch=std::min(s.batch_rays,n.area_rays.size());
    const auto cells=std::uint64_t(s.theta_count)*s.phi_count;
    if(batch>std::numeric_limits<unsigned>::max()*std::uint64_t(128))throw std::length_error("Atlas batch launch too large.");
    rays_.allocate(batch);audits_.allocate(batch);partial_.allocate(std::size_t(cells)*batch);
    mass_.allocate(cells);edges_.allocate(n.u_edges.size());
    rainbow::upload_host_to_device_sync(context_,edges_.address(),n.u_edges.data(),edges_.byte_size());
    pending_=true;mass_.zero_byte_async(context_.stream());
    Result r;
    for(std::size_t first=0;first<n.area_rays.size();first+=batch)
    {
        const auto count=std::min(batch,n.area_rays.size()-first);
        rainbow::upload_host_to_device_sync(context_,rays_.address(),n.area_rays.data()+first,count*sizeof(AreaRay));
        pending_=true;partial_.zero_byte_async(context_.stream());
        BatchParams p{};p.trace=n.trace;p.grid=n.view();p.grid.u_edges=edges_.data();
        p.rays=rays_.data();p.audits=audits_.data();p.mass=partial_.data();p.first=first;p.count=count;
        void* args[]={&p};launch("atlas_trace_area",args,static_cast<unsigned>((count+127)/128));sync();
        std::vector<RayAudit> audits(count);
        RAINBOW_CUDA_CHECK(cuMemcpyDtoH(audits.data(),audits_.address(),count*sizeof(RayAudit)));
        // A rejected batch is never committed to the accumulated table.
        add_audits(r,audits,{n.area_rays.data()+first,count},s.max_balance_error);
        const auto* in=partial_.data();auto* out=mass_.data();auto nr=std::uint64_t(count),nc=cells;
        void* reduce[]={&in,&out,&nc,&nr};launch("atlas_reduce_ordered",reduce,static_cast<unsigned>((cells+127)/128));
    }
    sync();std::vector<Sum> sums(cells);mass_.download(sums);r.mass_mm2.resize(cells);
    for(std::size_t k=0;k<cells;++k)r.mass_mm2[k]=sums[k].value;
    finish(r,s,n);return r;
}
}
