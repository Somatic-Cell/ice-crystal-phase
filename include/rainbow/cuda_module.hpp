#pragma once 

#include <rainbow/cuda_driver.hpp>

#include <cuda.h>

#include <filesystem>

namespace rainbow
{

class CudaModule final
{
public:
    explicit CudaModule(const CudaContext& cuda_context) noexcept;
    ~CudaModule() noexcept;

    CudaModule(const CudaModule&) = delete;
    CudaModule& operator=(const CudaModule&) = delete;
    CudaModule(CudaModule&&) = delete;
    CudaModule& operator=(CudaModule&&) = delete;

    // まだロードしていないときのみ使用可能
    void load_fatbin(const std::filesystem::path& path);

    // Kernel のシンボル名で関数を検索する
    // 見つからなければ例外を出す
    [[nodiscard]]
    CUfunction find_function(const char* name) const;

    void close();

    [[nodiscard]]
    bool close_noexcept() noexcept;

    [[nodiscard]]
    bool is_loaded() const noexcept {return module_ != nullptr; }

    [[nodiscard]]
    CUmodule handle() const noexcept {return module_; }

private:
    const CudaContext& cuda_context_;
    CUmodule module_ = nullptr;
    std::filesystem::path module_path_;
};
}