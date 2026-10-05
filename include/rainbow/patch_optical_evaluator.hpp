#pragma once

#include <rainbow/patch_optics_data.hpp>
#include <rainbow/focal_phase.hpp>
#include <rainbow/phase_interpolation.hpp>
#include <rainbow/polarization_transport.hpp>
#include <cmath>

namespace rainbow
{

// Per-hit optical arithmetic and per-direction accumulation, shared by CPU/CUDA.
// Stored Field32 includes Fresnel/interface phases, but not propagation phase.
// Parallel transport is our declared frame-identification convention, not an
// algorithm explicitly specified by Sadeghi et al. (2012).
struct PatchOpticalEvaluator
{
    struct Contribution
    {
        OpticalField64 path_field{};
        double incoherent_s = 0.0, incoherent_p = 0.0;
    };

    [[nodiscard]] HOST_DEVICE static std::uint32_t evaluate_patch(
        const OutgoingVertex* vertices, const std::uint32_t vertex_count,
        const OutgoingPatch& patch, const PatchQueryHit& hit,
        const Vec3 direction, const TransverseFrame target_frame,
        Contribution& output) noexcept
    {
        if(!vertices || (hit.flags & ~(BilinearBoundary | BilinearSingular | BilinearUsedFp64)) != 0u
           || hit.patch_id != patch.patch_id || hit.root_index > 1u
           || !(hit.u >= 0.0f && hit.u <= 1.0f && hit.v >= 0.0f && hit.v <= 1.0f)
           || !(hit.t > 0.0f && finite(hit.t)) || !(hit.residual >= 0.0f && finite(hit.residual)))
            return PatchOpticalInvalidInput;

        std::uint32_t pending = 0;
        if(patch.status == PatchCellStatus::NeedsRefinement) pending |= PatchOpticalRefinementPending;
        else if(!patch.has_regular_spherical_map()) return PatchOpticalInvalidInput;
        if(hit.flags & BilinearBoundary) pending |= PatchOpticalBoundaryPending;
        if(hit.flags & BilinearSingular) pending |= PatchOpticalSingularPending;
        if(pending != 0u) return pending; // no epsilon-area or zero-field substitution

        const double area = static_cast<double>(patch.incident_area_drop2);
        const double omega = ::fabs(static_cast<double>(patch.signed_solid_angle_sr));
        if(!(area > 0.0 && finite(area) && omega > 0.0 && finite(omega)))
            return PatchOpticalInvalidInput;

        Field32 transported[4]{};
        PhaseCycles paths[4]{};
        for(unsigned corner = 0; corner < 4; ++corner)
        {
            const auto index = patch.vertex_indices[corner];
            if(index >= vertex_count) return PatchOpticalInvalidInput;
            const auto& vertex = vertices[index];
            if(vertex.status != VertexStatus::Valid || !vertex.position_drop.is_finite()
               || !vertex.direction_drop.is_finite() || !vertex.basis_x.is_finite()
               || !finite_field(vertex.field)) return PatchOpticalInvalidInput;
            if(!vertex.optical_cycles.is_valid()) return PatchOpticalPhaseFailure;

            // Recover the stored right-handed transverse frame. Normalize only
            // the reconstructed second axis; the stored Jones components and
            // the first axis keep their existing representation.
            const auto second = vertex.direction_drop.cross(vertex.basis_x);
            const double length2 = second.cast<double>().dot(second.cast<double>());
            if(!(length2 > 0.0 && finite(length2))) return PatchOpticalFrameFailure;
            const TransverseFrame source_frame{vertex.basis_x, second.normalized()};
            if(!PolarizationTransport::try_transport(
                   vertex.direction_drop, source_frame, vertex.field,
                   direction, target_frame, transported[corner])) return PatchOpticalFrameFailure;
            paths[corner] = vertex.optical_cycles;
        }

        PhaseCycles interpolated_path{};
        if(!PhaseInterpolation::try_bilinear(paths[0], paths[1], paths[2], paths[3],
                                            hit.u, hit.v, interpolated_path))
            return PatchOpticalPhaseFailure;
        Complex32 phasor{};
        if(!interpolated_path.try_unit_phasor(phasor)) return PatchOpticalPhaseFailure;

        // Double interpolation/area scaling and accumulation; the existing
        // transport result and unit-phasor result are still FP32 by contract.
        const double amplitude = ::sqrt(area / omega);
        const double sr = amplitude * bilinear(transported[0].x.real, transported[1].x.real,
            transported[2].x.real, transported[3].x.real, hit.u, hit.v);
        const double si = amplitude * bilinear(transported[0].x.imag, transported[1].x.imag,
            transported[2].x.imag, transported[3].x.imag, hit.u, hit.v);
        const double pr = amplitude * bilinear(transported[0].y.real, transported[1].y.real,
            transported[2].y.real, transported[3].y.real, hit.u, hit.v);
        const double pi = amplitude * bilinear(transported[0].y.imag, transported[1].y.imag,
            transported[2].y.imag, transported[3].y.imag, hit.u, hit.v);
        const double c = static_cast<double>(phasor.real), s = static_cast<double>(phasor.imag);
        Contribution candidate{};
        candidate.incoherent_s = ::fma(sr, sr, si * si);
        candidate.incoherent_p = ::fma(pr, pr, pi * pi);
        candidate.path_field = {
            ::fma(sr, c, -si * s), ::fma(sr, s, si * c),
            ::fma(pr, c, -pi * s), ::fma(pr, s, pi * c)};
        if(!finite(candidate.incoherent_s) || !finite(candidate.incoherent_p)
           || !finite_field(candidate.path_field)) return PatchOpticalArithmeticFailure;
        output = candidate; // success-only update
        return PatchOpticalNone;
    }

    [[nodiscard]] HOST_DEVICE static PatchOpticalResult evaluate_direction(
        const PatchOpticsParams& input, const std::uint32_t index,
        const FocalPhaseConfig* focal_config = nullptr,
        FocalOpticalResult* focal_output = nullptr) noexcept
    {
        PatchOpticalResult result{};
        // Companion output is invalid until the entire known-hit evaluation finishes.
        // Early error exits must never leave a seemingly valid zero focal result.
        if(focal_output)
        {
            *focal_output = {};
            focal_output->flags = FocalInputError;
            invalidate_focal(*focal_output);
        }
        if(index >= input.direction_count || !input.directions || !input.offsets || !input.summaries)
        { result.flags = PatchOpticalInvalidInput; invalidate(result); return result; }

        const auto& summary = input.summaries[index];
        result.hit_count = summary.hit_count;
        result.query_flags = summary.flags;
        if(summary.flags & patch_query_error_mask)
        {
            result.flags = PatchOpticalQueryError;
            result.first_problem_patch_id = summary.first_problem_patch_id;
            result.rejected_hits = summary.hit_count;
            invalidate(result); return result;
        }
        const std::uint64_t begin = input.offsets[index], end = input.offsets[index + 1u];
        if(begin > end || end > input.hit_storage_count
           || std::uint64_t(summary.hit_count) > end - begin
           || std::uint64_t(summary.candidate_count) != end - begin
           || (summary.hit_count != 0u && (!input.hits || !input.patches || !input.vertices)))
        { result.flags = PatchOpticalInvalidInput; result.rejected_hits = summary.hit_count;
          invalidate(result); return result; }

        // Pending flags in the summary are also preserved, even if supplied
        // externally without a corresponding record. Never certify such input.
        if(summary.flags & PatchQueryRefinementHit) result.flags |= PatchOpticalRefinementPending;
        if(summary.flags & PatchQueryBoundaryHit) result.flags |= PatchOpticalBoundaryPending;
        if(summary.flags & PatchQuerySingularHit) result.flags |= PatchOpticalSingularPending;
        TransverseFrame target_frame{};
        if(!PolarizationTransport::try_make_scattering_frame(input.incident_direction,
               input.incident_basis_x, input.directions[index], target_frame))
        { result.flags |= PatchOpticalFrameFailure; result.rejected_hits = summary.hit_count;
          invalidate(result); return result; }

        Sum sr{}, si{}, pr{}, pi{}, incoherent_s{}, incoherent_p{};
        const bool with_focal = focal_config && focal_output;
        FocalOpticalResult focal{};
        Sum focal_sr{}, focal_si{}, focal_pr{}, focal_pi{}, family_envelope[4]{};
        if(with_focal && !FocalPhase::layout_valid(*focal_config)) focal.flags |= FocalInvalidGeometry;

        for(std::uint32_t i = 0; i < summary.hit_count; ++i)
        {
            const auto& hit = input.hits[begin + i];
            std::uint32_t flags = 0;
            Contribution contribution{};
            // The upstream query sorts and removes exact duplicate reports.
            // Reject a broken ordering/duplicate contract instead of summing twice.
            if(i != 0u && !strictly_before(input.hits[begin + i - 1u], hit))
                flags = PatchOpticalInvalidInput;
            else if(hit.compact_index >= input.patch_count) flags = PatchOpticalInvalidInput;
            else flags = evaluate_patch(input.vertices, input.vertex_count,
                    input.patches[hit.compact_index], hit, input.directions[index], target_frame, contribution);

            if(flags != 0u)
            {
                result.flags |= flags; ++result.rejected_hits;
                result.refinement_hits += (flags & PatchOpticalRefinementPending) != 0u;
                result.boundary_hits += (flags & PatchOpticalBoundaryPending) != 0u;
                result.singular_hits += (flags & PatchOpticalSingularPending) != 0u;
                if(hit.patch_id < result.first_problem_patch_id) result.first_problem_patch_id = hit.patch_id;
                continue;
            }
            ++result.evaluated_hits;
            sr.add(contribution.path_field.s_real); si.add(contribution.path_field.s_imag);
            pr.add(contribution.path_field.p_real); pi.add(contribution.path_field.p_imag);
            incoherent_s.add(contribution.incoherent_s); incoherent_p.add(contribution.incoherent_p);
            if(with_focal)
            {
                FocalPhaseEstimate estimate{};
                const auto code = FocalPhase::estimate(input.vertices, input.vertex_count,
                    input.patches[hit.compact_index], input.incident_direction, *focal_config, estimate);
                if(code != FocalNone)
                {
                    focal.flags |= code;
                    if(hit.patch_id < focal.first_problem_patch_id) focal.first_problem_patch_id = hit.patch_id;
                }
                else
                {
                    const auto f = FocalPhase::apply(contribution.path_field, estimate.quarter_turns);
                    focal_sr.add(f.s_real); focal_si.add(f.s_imag);
                    focal_pr.add(f.p_real); focal_pi.add(f.p_imag);
                    family_envelope[estimate.family].add(contribution.incoherent_s + contribution.incoherent_p);
                    ++focal.family_hits[estimate.family];
                    ++focal.corrected_hits;
                    focal.extra_quarter_turn_hits += estimate.extra_quarter_turn;
                }
            }

        }
        result.regular_partial_path_field = {sr.value(), si.value(), pr.value(), pi.value()};
        result.regular_partial_incoherent_s = incoherent_s.value();
        result.regular_partial_incoherent_p = incoherent_p.value();
        result.regular_partial_path_s = ::fma(sr.value(), sr.value(), si.value() * si.value());
        result.regular_partial_path_p = ::fma(pr.value(), pr.value(), pi.value() * pi.value());
        if(!finite_field(result.regular_partial_path_field)
           || !finite(result.regular_partial_incoherent_s) || !finite(result.regular_partial_incoherent_p)
           || !finite(result.regular_partial_path_s) || !finite(result.regular_partial_path_p))
            result.flags |= PatchOpticalArithmeticFailure;
        if(result.flags & patch_optical_error_mask) invalidate(result);
        if(with_focal)
        {
            if(result.flags & patch_optical_error_mask) focal.flags |= FocalInputError;
            if(result.flags & patch_optical_pending_mask) focal.flags |= FocalInputPending;
            if(result.first_problem_patch_id < focal.first_problem_patch_id)
                focal.first_problem_patch_id = result.first_problem_patch_id;
            focal.field = {focal_sr.value(), focal_si.value(), focal_pr.value(), focal_pi.value()};
            focal.intensity_s = ::fma(focal.field.s_real, focal.field.s_real, focal.field.s_imag * focal.field.s_imag);
            focal.intensity_p = ::fma(focal.field.p_real, focal.field.p_real, focal.field.p_imag * focal.field.p_imag);
            for(unsigned f = 0; f < 4; ++f) focal.family_incoherent[f] = family_envelope[f].value();
            if(!finite_field(focal.field) || !finite(focal.intensity_s + focal.intensity_p))
                focal.flags |= FocalArithmeticError;
            // No complete-looking focal intensity assembled from only some hits.
            if(!focal.valid()) invalidate_focal(focal);
            *focal_output = focal;
        }
        return result;
    }

private:
    HOST_DEVICE static void invalidate_focal(FocalOpticalResult& r) noexcept
    {
        r.field = {patch_optical_nan, patch_optical_nan, patch_optical_nan, patch_optical_nan};
        r.intensity_s = r.intensity_p = patch_optical_nan;
    }
    // Neumaier compensated summation. No atomic operations; one thread owns
    // one direction and consumes every hit in the existing deterministic order.
    struct Sum
    {
        double sum = 0.0, correction = 0.0;
        HOST_DEVICE void add(const double x) noexcept
        {
            const double next = sum + x;
            correction += ::fabs(sum) >= ::fabs(x) ? (sum - next) + x : (x - next) + sum;
            sum = next;
        }
        [[nodiscard]] HOST_DEVICE double value() const noexcept { return sum + correction; }
    };
    [[nodiscard]] HOST_DEVICE static double bilinear(
        const double a, const double b, const double c, const double d,
        const double u, const double v) noexcept
    {
        const double lower = ::fma(u, b - a, a), upper = ::fma(u, d - c, c);
        return ::fma(v, upper - lower, lower);
    }
    [[nodiscard]] HOST_DEVICE static bool finite(const double x) noexcept
    { return x >= -0x1.fffffffffffffp1023 && x <= 0x1.fffffffffffffp1023; }
    [[nodiscard]] HOST_DEVICE static bool finite_field(const Field32 f) noexcept
    { return finite(f.x.real) && finite(f.x.imag) && finite(f.y.real) && finite(f.y.imag); }
    [[nodiscard]] HOST_DEVICE static bool finite_field(const OpticalField64 f) noexcept
    { return finite(f.s_real) && finite(f.s_imag) && finite(f.p_real) && finite(f.p_imag); }
    [[nodiscard]] HOST_DEVICE static bool strictly_before(const PatchQueryHit& a, const PatchQueryHit& b) noexcept
    {
        if(a.patch_id != b.patch_id) return a.patch_id < b.patch_id;
        if(a.compact_index != b.compact_index) return a.compact_index < b.compact_index;
        return a.root_index < b.root_index;
    }
    HOST_DEVICE static void invalidate(PatchOpticalResult& result) noexcept
    {
        result.regular_partial_path_field = {patch_optical_nan, patch_optical_nan,
                                             patch_optical_nan, patch_optical_nan};
        result.regular_partial_incoherent_s = result.regular_partial_incoherent_p = patch_optical_nan;
        result.regular_partial_path_s = result.regular_partial_path_p = patch_optical_nan;
    }
};

} // namespace rainbow
