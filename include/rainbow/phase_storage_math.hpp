#pragma once

#include <cmath>
#include <cstdint>

#if defined(__CUDACC__)
#define RAINBOW_STORAGE_HD __host__ __device__
#else
#define RAINBOW_STORAGE_HD
#endif

namespace rainbow
{
// Numerical kernels shared by CPU regression tests and the CUDA entry points.
// No CUDA Runtime API, random generation, or neural-network dependency.
namespace phase_storage_math
{
inline constexpr double pi = 3.141592653589793238462643383279502884;
struct Sum
{
    double s = 0, c = 0;
    RAINBOW_STORAGE_HD void add(double x) noexcept
    {
        const double t = s + x;
        c += ::fabs(s) >= ::fabs(x) ? (s-t)+x : (x-t)+s;
        s = t;
    }
    [[nodiscard]] RAINBOW_STORAGE_HD double value() const noexcept { return s+c; }
};
struct Grid
{
    const double* values = nullptr; // theta-major density, NOT CDF or log intensity
    const double* edges = nullptr;  // u=(1-cos(theta))/2, length nt+1
    std::uint32_t nt = 0, np = 0;
};
[[nodiscard]] RAINBOW_STORAGE_HD inline double mean_cosine(double lower_u, double upper_u) noexcept
{
    // Exact cell average for density constant with respect to solid angle.
    return 1.0 - (lower_u + upper_u);
}
[[nodiscard]] RAINBOW_STORAGE_HD inline int wrap(int i, int n) noexcept
{
    const int r = i % n; return r < 0 ? r+n : r;
}
[[nodiscard]] RAINBOW_STORAGE_HD inline double angular_distance(double a, double b, double dp) noexcept
{
    const double x = ::sin(0.5*(a-b)), y = ::sin(0.5*dp);
    const double h = ::fma(::sin(a)*::sin(b), y*y, x*x);
    return 2.0*::asin(::sqrt(::fmin(1.0, ::fmax(0.0, h))));
}
// Direct, solid-angle-weighted spherical Gaussian quadrature. Not a separable
// image-space Gaussian. Longitude is periodic; polar caps may include all phi.
// Each output row is normalized: constants are preserved. On a finite angular
// grid this does NOT guarantee exact global mass conservation. The generator
// records that integral change and normalizes only when constructing the CDF.
[[nodiscard]] RAINBOW_STORAGE_HD inline double gaussian_at(
    Grid g, std::uint32_t row, std::uint32_t col, double sigma, double support) noexcept
{
    const double dt = pi/double(g.nt), dp = 2*pi/double(g.np);
    const double t = (double(row)+0.5)*dt, radius = support*sigma;
    int low = int(::floor((t-radius)/dt))-1;
    int high = int(::ceil((t+radius)/dt))+1;
    low = low < 0 ? 0 : low; high = high >= int(g.nt) ? int(g.nt)-1 : high;
    Sum weighted{}, weights{};
    for(int r = low; r <= high; ++r)
    {
        const double tr = (double(r)+0.5)*dt;
        if(::fabs(t-tr) > radius+0x1p-44) continue;
        const double st = ::sin(t)*::sin(tr);
        const double x = st > 0 ? (::cos(radius)-::cos(t)*::cos(tr))/st : -1.0;
        int width;
        if(x <= -1.0) width = int(g.np);
        else if(x >= 1.0) width = 1;
        else width = int(::ceil(::acos(x)/dp))+1;
        const int count = 2*width+1 < int(g.np) ? 2*width+1 : int(g.np);
        const int start = count == int(g.np) ? 0 : int(col)-width;
        const double du = g.edges[r+1]-g.edges[r];
        for(int j = 0; j < count; ++j)
        {
            const int c = wrap(start+j, int(g.np));
            const double distance = angular_distance(t, tr, double(c-int(col))*dp);
            if(distance > radius+0x1p-44) continue;
            const double s = distance/sigma;
            const double w = ::exp(-0.5*s*s)*du; // common dphi cancels
            weights.add(w);
            weighted.add(w*g.values[std::uint64_t(r)*g.np+std::uint32_t(c)]);
        }
    }
    return weighted.value()/weights.value(); // center guarantees positive weight
}
// Nested grids only. Preserve the mass of the supplied piecewise-constant
// density; theta interpolation/averaging without its solid angle is incorrect.
[[nodiscard]] RAINBOW_STORAGE_HD inline double coarsen_at(
    Grid g, std::uint32_t row, std::uint32_t col,
    std::uint32_t out_nt, std::uint32_t out_np) noexcept
{
    const auto rt = g.nt/out_nt, rp = g.np/out_np;
    const auto begin = row*rt, end = (row+1u)*rt;
    Sum weighted{};
    for(auto r=begin; r<end; ++r)
        for(auto c=col*rp; c<(col+1u)*rp; ++c)
            weighted.add(g.values[std::uint64_t(r)*g.np+c]*(g.edges[r+1]-g.edges[r]));
    return weighted.value()/((g.edges[end]-g.edges[begin])*double(rp));
}
} // namespace phase_storage_math
} // namespace rainbow
#undef RAINBOW_STORAGE_HD
