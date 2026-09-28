#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#endif

#include <rainbow/optix_error.hpp>


#include <optix_stubs.h>

#include <cstdio>
#include <sstream>
#include <stdexcept>

namespace rainbow::detail
{

void check_optix(
    const OptixResult result,
    const char* expression,
    const char* file,
    const int line)
{
    if(result == OPTIX_SUCCESS)
    {
        return;
    }

    // optix_stubs.h は，optixInit() の失敗時にもエラー名を取得できる実装になっている
    const char* error_name = optixGetErrorName(result);
    const char* error_description = optixGetErrorString(result);
    
    std::ostringstream message;
    message
        << expression
        << " failed at "
        << file
        << ':'
        << line
        << '\n'
        << "OptiX error: "
        << (error_name != nullptr ? error_name : "unknown")
        << " (" << static_cast<int>(result) << ")"
        << '\n'
        << "Description: "
        << (error_description != nullptr ? error_description : "unavailable");
    throw std::runtime_error(message.str());
}

bool report_optix_cleanup_result(
    const OptixResult result,
    const char* expression) noexcept
{
    if(result == OPTIX_SUCCESS)
    {
        return true;
    }

    const char* error_name = optixGetErrorName(result);
    std::fprintf(
        stderr, "[OptiX cleanup] %s: %s: (%d)\n",
        expression,
        error_name != nullptr ? error_name : "unknown",
        static_cast<int>(result)
    );
    return false;
}

}