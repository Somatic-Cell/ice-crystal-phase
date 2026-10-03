#pragma once

#include <rainbow/field32.hpp>
#include <rainbow/vec3.hpp>

#include <cmath>

namespace rainbow
{

// A documented convention for comparing Jones vectors carried by neighbouring
// rays whose propagation directions differ.
//
// The paper states that each ray carries a transverse frame but does not specify
// how four different frames are reconciled during patch interpolation.  Here a
// source field is parallel-transported to the target direction using the unique
// minimum rotation that maps source_direction to target_direction.  The
// transported physical vector is then expressed in one common target frame.
// This makes the result invariant to a mere rotation (gauge change) of either
// ray's stored transverse basis.
struct PolarizationTransport
{
    [[nodiscard]]
    HOST_DEVICE static bool try_make_scattering_frame(
        const Vec3 incident_direction,
        const Vec3 incident_basis_x,
        const Vec3 outgoing_direction,
        TransverseFrame& output) noexcept
    {
        if(!is_unit(incident_direction)
           || !is_unit(outgoing_direction)
           || !is_unit(incident_basis_x)
           || ::fabsf(incident_direction.dot(incident_basis_x)) > frame_tolerance)
        {
            return false;
        }

        Vec3 e0 = incident_direction.cross(outgoing_direction);
        float e0_norm2 = e0.dot(e0);

        if(e0_norm2 > degeneracy_threshold)
        {
            e0 = e0 / ::sqrtf(e0_norm2);
        }
        else
        {
            // At exact forward/backward scattering the scattering plane is
            // undefined.  Use the transported incident reference axis only as
            // a deterministic gauge; QueryDirectionGrid avoids the exact poles.
            e0 = incident_basis_x
                - outgoing_direction * outgoing_direction.dot(incident_basis_x);
            e0_norm2 = e0.dot(e0);
            if(!(e0_norm2 > degeneracy_threshold))
            {
                return false;
            }
            e0 = e0 / ::sqrtf(e0_norm2);
        }

        Vec3 e1 = outgoing_direction.cross(e0);
        const float e1_norm2 = e1.dot(e1);
        if(!(e1_norm2 > degeneracy_threshold))
        {
            return false;
        }
        e1 = e1 / ::sqrtf(e1_norm2);

        const TransverseFrame candidate{e0, e1};
        if(!is_frame(outgoing_direction, candidate))
        {
            return false;
        }
        output = candidate;
        return true;
    }

    [[nodiscard]]
    HOST_DEVICE static bool try_transport(
        const Vec3 source_direction,
        const TransverseFrame source_frame,
        const Field32 source_field,
        const Vec3 target_direction,
        const TransverseFrame target_frame,
        Field32& output) noexcept
    {
        if(!is_unit(source_direction)
           || !is_unit(target_direction)
           || !is_frame(source_direction, source_frame)
           || !is_frame(target_direction, target_frame))
        {
            return false;
        }

        // FP32 の単位長誤差を持ち越さないよう，変換後に FP64 で再正規化する．
        // 最小回転の式と反対方向の判定は，単位ベクトルを前提にする．
        const DVec3 source_w = source_direction.cast<double>().normalized();
        const DVec3 target_w = target_direction.cast<double>().normalized();
        const double cosine = source_w.dot(target_w);

        // There is no unique shortest rotation for antipodal directions.
        // A regular patch should never ask us to transport across such a span.
        if(!(1.0 + cosine > antipodal_threshold))
        {
            return false;
        }

        const DVec3 source_e0 = source_frame.e0.cast<double>();
        const DVec3 source_e1 = source_frame.e1.cast<double>();
        const DVec3 real_world =
            source_e0 * static_cast<double>(source_field.x.real)
            + source_e1 * static_cast<double>(source_field.y.real);
        const DVec3 imag_world =
            source_e0 * static_cast<double>(source_field.x.imag)
            + source_e1 * static_cast<double>(source_field.y.imag);

        const DVec3 real_transported =
            minimum_rotation(real_world, source_w, target_w, cosine);
        const DVec3 imag_transported =
            minimum_rotation(imag_world, source_w, target_w, cosine);
        const DVec3 target_e0 = target_frame.e0.cast<double>();
        const DVec3 target_e1 = target_frame.e1.cast<double>();

        const double x_real = real_transported.dot(target_e0);
        const double x_imag = imag_transported.dot(target_e0);
        const double y_real = real_transported.dot(target_e1);
        const double y_imag = imag_transported.dot(target_e1);
        if(!fits_float(x_real) || !fits_float(x_imag)
           || !fits_float(y_real) || !fits_float(y_imag))
        {
            return false;
        }

        output = {
            {static_cast<float>(x_real), static_cast<float>(x_imag)},
            {static_cast<float>(y_real), static_cast<float>(y_imag)}
        };
        return true;
    }

private:
    using DVec3 = Vec3T<double>;

    static constexpr float frame_tolerance = 0x1p-10f;
    static constexpr float degeneracy_threshold = 0x1p-30f;
    static constexpr double antipodal_threshold = 0x1p-30;

    [[nodiscard]]
    HOST_DEVICE static bool is_unit(const Vec3 value) noexcept
    {
        if(!value.is_finite()) return false;
        return ::fabsf(value.dot(value) - 1.0f) <= frame_tolerance;
    }

    [[nodiscard]]
    HOST_DEVICE static bool is_frame(
        const Vec3 direction,
        const TransverseFrame frame) noexcept
    {
        if(!is_unit(frame.e0) || !is_unit(frame.e1)) return false;
        if(::fabsf(direction.dot(frame.e0)) > frame_tolerance
           || ::fabsf(direction.dot(frame.e1)) > frame_tolerance
           || ::fabsf(frame.e0.dot(frame.e1)) > frame_tolerance)
        {
            return false;
        }
        return frame.e0.cross(frame.e1).dot(direction) >= 1.0f - 4.0f * frame_tolerance;
    }

    [[nodiscard]]
    HOST_DEVICE static DVec3 minimum_rotation(
        const DVec3 value,
        const DVec3 source_direction,
        const DVec3 target_direction,
        const double cosine) noexcept
    {
        // Rodrigues without evaluating an angle:
        // R x = x + v x x + v x (v x x)/(1 + a.b), v = a x b.
        const DVec3 v = source_direction.cross(target_direction);
        const DVec3 vx = v.cross(value);
        return value + vx + v.cross(vx) / (1.0 + cosine);
    }

    [[nodiscard]]
    HOST_DEVICE static bool fits_float(const double value) noexcept
    {
        constexpr double maximum = 0x1.fffffep127;
        return value >= -maximum && value <= maximum;
    }
};

} // namespace rainbow
