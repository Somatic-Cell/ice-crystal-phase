#pragma once

#include <optix.h>

namespace rainbow::detail
{

void check_optix(
    const OptixResult result,
    const char* expression,
    const char* file,
    const int line);
 
[[nodiscard]]
bool report_optix_cleanup_result(
    const OptixResult result,
    const char* expression) noexcept;

} // namespace rainbow::detail

#define RAINBOW_OPTIX_CHECK(expression) \
    do \
    { \
        ::rainbow::detail::check_optix( \
            (expression), #expression, __FILE__, __LINE__ \
        ); \
    } while (false)
    