#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include <windows.h>

#include <rainbow/cuda_driver.hpp>
#include <rainbow/optix.hpp>

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>

namespace{
[[nodiscard]]
std::filesystem::path executable_directory()
{
    std::wstring path_buffer(256, L'\0');

    for(;;)
    {
        const DWORD character_count = GetModuleFileNameW(
            nullptr,
            path_buffer.data(),
            static_cast<DWORD>(path_buffer.size())
        );

        if(character_count == 0)
        {
            const DWORD error_code = GetLastError();
            throw std::system_error(
                static_cast<int>(error_code),
                std::system_category(),
                "GetModuleFileNameW failed");
        }

        if(character_count < path_buffer.size())
        {
            path_buffer.resize(character_count);
            return std::filesystem::path{path_buffer}.parent_path();
        }

        // 次回の DWORD への変換と，size の倍増がオーバーフローしないようにする．
        if(path_buffer.size() > (std::numeric_limits<DWORD>::max)() / 2)
        {
            throw std::length_error("Executable path buffer is too large.");
        }
        path_buffer.resize(path_buffer.size() * 2);

    }
}

} // namespace

int main(
    const int argc,
    char* argv[])
{
    try
    {
        // 従来どおり，省略可能な引数は fatbin と OptiX IR の順に受け取る．
        if(argc < 1 || argc > 3)
        {
            std::cerr
                << "Usage: rainbow.exe [fatbin_path [optixir_path]]\n";
            return EXIT_FAILURE;
        }

        const std::filesystem::path modules_directory =
            executable_directory() / "modules";

        // 引数なし: .exe の隣の modules/ を使う．
        // 引数あり: 指定されたパスを使う．相対パスは起動時の作業ディレクトリ基準．
        const std::filesystem::path fatbin_path =
            argc >= 2
                ? std::filesystem::absolute(std::filesystem::path{argv[1]})
                : modules_directory / "patch_build.fatbin";

        const std::filesystem::path optixir_path =
            argc >= 3
                ? std::filesystem::absolute(std::filesystem::path{argv[2]})
                : modules_directory / "optix_smoke.optixir";

        // 所有者は一つだけ作り，M1 と M2 の両方より長く生存させる．
        // ordinal 0 は CUDA から見える最初の device を指定する．
        rainbow::CudaContext cuda_context{0};

        // 両関数は同じ context / stream を参照で借りる．コピーも移動もしない．
        // M1 が例外で失敗した場合は，M2 へ進まず catch へ移る．
        rainbow::run_cuda_driver_smoke_test(cuda_context, fatbin_path);
        rainbow::run_optix_smoke_test(cuda_context, optixir_path);

        // 関数内の一時リソースは各テストが解放済み．このスコープを抜けると
        // cuda_context が stream と primary context の retain 参照を解放する．
        return EXIT_SUCCESS;
    }
    catch(const std::exception& exception)
    {
        std::cerr
            << "Error: "
            << exception.what()
            << '\n';

        return 1;
    }
    catch(...)
    {
        std::cerr << "Error: an unknown exception was thrown.\n";
        return EXIT_FAILURE;
    }
}