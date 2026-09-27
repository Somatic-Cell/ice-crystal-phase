#pragma once

#include <cstdint>

namespace rainbow
{
struct alignas(8) TraceLaunchParams{
    std::uint32_t* output;
    std::uint32_t count;
    std::uint32_t reserved;
};

static_assert(sizeof(TraceLaunchParams) == 16);

static_assert(alignof(TraceLaunchParams) == 8);

}