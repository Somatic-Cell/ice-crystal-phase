#pragma once

#include <rainbow/raindrop_trace_data.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace rainbow
{
// theta は入射光の「進行方向」と出射方向のなす角．phi=0 は incident_basis_x．
// セル中心を使い，theta=0,pi と phi の seam を重複サンプリングしない．
// 非等面積の longitude/latitude 格子．積分時は solid_angle() の重みが必要．
struct QueryDirectionGrid
{
    std::uint32_t theta_count=0, phi_count=0;
    Vec3 axis{}, e0{}, e1{};

    [[nodiscard]] static QueryDirectionGrid make(
        const RaindropTraceConfig& config, const std::uint32_t ntheta,
        const std::uint32_t nphi)
    {
        if(ntheta==0 || nphi==0 || std::uint64_t(ntheta)*nphi>0xffffffffull)
            throw std::invalid_argument("Invalid query grid dimensions.");
        const auto w=config.incident_direction;
        const auto ex=config.incident_basis_x;
        if(!w.is_finite() || !ex.is_finite()
           || std::abs(w.dot(w)-1.0f)>1e-5f || std::abs(ex.dot(ex)-1.0f)>1e-5f
           || std::abs(w.dot(ex))>1e-5f)
            throw std::invalid_argument("Query grid requires an orthonormal incident frame.");
        return {ntheta,nphi,w,ex,w.cross(ex).normalized()};
    }
    [[nodiscard]] std::uint64_t size() const noexcept {return std::uint64_t(theta_count)*phi_count;}
    [[nodiscard]] double theta(const std::uint32_t row) const noexcept
    {return pi*(double(row)+0.5)/double(theta_count);}
    [[nodiscard]] double phi(const std::uint32_t column) const noexcept
    {return -pi+2.0*pi*(double(column)+0.5)/double(phi_count);}
    [[nodiscard]] double solid_angle(const std::uint32_t row) const noexcept
    {
        const double a=pi*double(row)/double(theta_count);
        const double b=pi*double(row+1)/double(theta_count);
        return (2.0*pi/double(phi_count))*(std::cos(a)-std::cos(b));
    }
    [[nodiscard]] std::vector<Vec3> directions() const
    {
        std::vector<Vec3> output(static_cast<std::size_t>(size()));
        for(std::uint32_t row=0;row<theta_count;++row)
        {
            const double t=theta(row),st=std::sin(t),ct=std::cos(t);
            for(std::uint32_t col=0;col<phi_count;++col)
            {
                const double p=phi(col),x=st*std::cos(p),y=st*std::sin(p);
                const Vec3 w{
                    float(double(axis.x)*ct+double(e0.x)*x+double(e1.x)*y),
                    float(double(axis.y)*ct+double(e0.y)*x+double(e1.y)*y),
                    float(double(axis.z)*ct+double(e0.z)*x+double(e1.z)*y)};
                output[std::size_t(row)*phi_count+col]=w.normalized();
            }
        }
        return output;
    }
    static constexpr double pi=3.141592653589793238462643383279502884;
};
} // namespace rainbow
