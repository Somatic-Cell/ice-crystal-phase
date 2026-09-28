#include <rainbow/cuda_module.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/read_binary_file.hpp>

#include <iostream>
#include <stdexcept>
#include <vector>

namespace rainbow
{

CudaModule::CudaModule(const CudaContext& cuda_context) noexcept
    : cuda_context_(cuda_context)
{
}

CudaModule::~CudaModule() noexcept
{
    static_cast<void>(close_noexcept());
}

void CudaModule::load_fatbin(const std::filesystem::path& path)
{
    if(is_loaded())
    {
        throw std::logic_error("CudaModule is already loaded.");
    }

    // path の複製と読み込みは例外を出し得るため，GPU module 取得より先に行う．
    std::filesystem::path new_path = path;
    const std::vector<char> fatbin = read_binary_file(path);
    cuda_context_.make_current();

    CUmodule new_module = nullptr;
    RAINBOW_CUDA_CHECK(cuModuleLoadData(&new_module, fatbin.data()));

    // 成功後の状態更新では，例外を出さない swap と raw handle の代入だけを使う．
    module_path_.swap(new_path);
    module_ = new_module;
}

CUfunction CudaModule::find_function(const char* name) const
{
    if(!is_loaded())
    {
        throw std::logic_error("Cannot look up a function in an unloaded CudaModule.");
    }
    if(name == nullptr || name[0] == '\0')
    {
        throw std::invalid_argument("CUDA kernel name must not be empty.");
    }

    cuda_context_.make_current();
    CUfunction function = nullptr;
    const CUresult result = cuModuleGetFunction(&function, module_, name);
    if(result == CUDA_ERROR_NOT_FOUND)
    {
        // 以前の診断を維持する．モジュール配置先の更新漏れも切り分けられる．
        std::cerr << "CUDA module: " << module_path_
                  << "\nExpected kernel: " << name
                  << "\nCheck extern \"C\" and rebuild / update the staged .fatbin file.\n";
    }
    detail::check_cuda(result, "cuModuleGetFunction", __FILE__, __LINE__);
    return function;
}

void CudaModule::close()
{
    if(!close_noexcept())
    {
        throw std::runtime_error("CudaModule cleanup failed. See stderr.");
    }
}

bool CudaModule::close_noexcept() noexcept
{
    if(!is_loaded())
    {
        return true;
    }
    if(!detail::report_cuda_cleanup_result(
           cuCtxSetCurrent(cuda_context_.handle()), "cuCtxSetCurrent(module cleanup)"))
    {
        return false;
    }

    const CUresult result = cuModuleUnload(module_);
    module_ = nullptr;
    module_path_.clear();
    return detail::report_cuda_cleanup_result(result, "cuModuleUnload(CudaModule)");
}

} // namespace rainbow
