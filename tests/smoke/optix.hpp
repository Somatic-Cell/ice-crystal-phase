#pragma once

#include <filesystem>

namespace rainbow
{
    class CudaContext;

    void run_optix_smoke_test(
        const CudaContext& cuda_context,
        const std::filesystem::path& optixir_path
    );
}