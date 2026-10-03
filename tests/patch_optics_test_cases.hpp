#pragma once

#include <rainbow/patch_optical_evaluator.hpp>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace rainbow::tests
{
inline void require_optics(const bool ok, const char* message)
{ if(!ok) throw std::runtime_error(message); }
inline bool near_optics(const double a, const double b, const double scale = 1.0)
{
    return std::isfinite(a) && std::isfinite(b)
        && std::abs(a - b) <= 8e-6 * (scale + std::abs(b));
}

struct OpticalExpectation
{
    double incoherent_s = 0, incoherent_p = 0, path_s = 0, path_p = 0;
    OpticalField64 field{};
    std::uint32_t flags = 0, evaluated = 0, rejected = 0;
};

// Arithmetic fixtures, deliberately independent of GAS construction. Their
// assigned solid angles test the evaluator, not the geometry's area estimator.
// Integration with real curved patches is tested separately on GPU.
struct OpticsFixture
{
    std::vector<OutgoingVertex> vertices;
    std::vector<OutgoingPatch> patches;
    std::vector<Vec3> directions;
    std::vector<std::uint64_t> offsets{0};
    std::vector<PatchQuerySummary> summaries;
    std::vector<PatchQueryHit> hits;
    std::vector<OpticalExpectation> expected;
    Vec3 incident{1,0,0}, incident_basis{0,0,1};

    [[nodiscard]] PatchOpticsParams view() const noexcept
    {
        PatchOpticsParams p{};
        p.vertices=vertices.data();p.patches=patches.data();p.directions=directions.data();
        p.offsets=offsets.data();p.summaries=summaries.data();p.hits=hits.data();
        p.hit_storage_count=hits.size();p.vertex_count=static_cast<std::uint32_t>(vertices.size());
        p.patch_count=static_cast<std::uint32_t>(patches.size());
        p.direction_count=static_cast<std::uint32_t>(directions.size());
        p.incident_direction=incident;p.incident_basis_x=incident_basis;
        return p;
    }

    void add_hit(const Field32 field, const PhaseCycles phase, const float area = 1.0f,
                 const float omega = 1.0f, const std::uint32_t hit_flags = 0,
                 const PatchCellStatus status = PatchCellStatus::RegularPositive)
    {
        const auto pi=static_cast<std::uint32_t>(patches.size());
        OutgoingPatch patch{};patch.patch_id=pi;patch.status=status;
        patch.incident_area_drop2=area;patch.signed_solid_angle_sr=omega;
        for(unsigned c=0;c<4;++c)
        {
            patch.vertex_indices[c]=static_cast<std::uint32_t>(vertices.size());
            OutgoingVertex v{};v.direction_drop={0,0,1};v.basis_x={0,-1,0};
            v.field=field;v.optical_cycles=phase;v.status=VertexStatus::Valid;
            vertices.push_back(v);
        }
        patches.push_back(patch);
        hits.push_back({pi,pi,0,hit_flags,.5f,.5f,1.0f,0.0f});
    }
    void end_direction(OpticalExpectation e, const std::uint32_t query_flags = 0,
                       const std::uint32_t unused_slots = 0)
    {
        const auto count=static_cast<std::uint32_t>(hits.size()-offsets.back());
        directions.push_back({0,0,1});
        PatchQuerySummary summary{};summary.hit_count=count;
        summary.candidate_count=count+unused_slots;summary.flags=query_flags;
        summaries.push_back(summary);expected.push_back(e);
        // Trailing capacity contains invalid records. A correct evaluator never
        // reads them; there is no fixed upper bound on the real hit list.
        for(std::uint32_t i=0;i<unused_slots;++i)
            hits.push_back({0xffffffffu,0xffffffffu,0,0,-1,-1,-1,-1});
        offsets.push_back(hits.size());
    }
};

inline OpticsFixture make_optics_fixture()
{
    OpticsFixture f;
    const Field32 one{{1,0},{0,0}};
    f.add_hit(one,{0,0}); f.end_direction({1,0,1,0,{1,0,0,0},0,1,0});
    f.add_hit(one,{0,0}); f.add_hit(one,{1000,0});
    f.end_direction({2,0,4,0,{2,0,0,0},0,2,0},0,3);
    f.add_hit(one,{0,0}); f.add_hit(one,{1000,.5f});
    f.end_direction({2,0,0,0,{0,0,0,0},0,2,0});
    f.add_hit(one,{0,0}); f.add_hit(one,{1000,.25f});
    f.end_direction({2,0,2,0,{1,1,0,0},0,2,0});
    f.add_hit(one,{0,0}); f.add_hit({{0,0},{1,0}},{0,.5f});
    f.end_direction({1,1,1,1,{1,0,-1,0},0,2,0});
    f.add_hit(one,{0,.75f},2,-.5f,0,PatchCellStatus::RegularNegative);
    f.end_direction({4,0,4,0,{0,-2,0,0},0,1,0});
    f.end_direction({}); // actual zero intensity, not an invalid direction
    f.add_hit({},{0,0});f.end_direction({0,0,0,0,{},0,1,0});
    // Compensation preserves the small wave between cancelling large waves.
    const float large=1e20f;const double ld=large;
    f.add_hit({{large,0},{0,0}},{0,0});f.add_hit(one,{0,0});
    f.add_hit({{-large,0},{0,0}},{0,0});
    f.end_direction({2*ld*ld+1,0,1,0,{1,0,0,0},0,3,0});
    f.add_hit(one,{0,0});f.add_hit(one,{0,0},1,0,0,PatchCellStatus::NeedsRefinement);
    f.end_direction({1,0,1,0,{1,0,0,0},PatchOpticalRefinementPending,1,1},PatchQueryRefinementHit);
    f.add_hit(one,{0,0});f.add_hit(one,{0,0},1,1,BilinearBoundary);
    f.end_direction({1,0,1,0,{1,0,0,0},PatchOpticalBoundaryPending,1,1},PatchQueryBoundaryHit);
    f.add_hit(one,{0,0},1,1,BilinearSingular);
    f.end_direction({0,0,0,0,{},PatchOpticalSingularPending,0,1},PatchQuerySingularHit);
    f.add_hit(one,{0,0},1,0);
    f.end_direction({0,0,0,0,{},PatchOpticalInvalidInput,0,1});
    f.add_hit(one,{0,1});
    f.end_direction({0,0,0,0,{},PatchOpticalPhaseFailure,0,1});
    f.add_hit(one,{0,0});f.hits.back().u=-.01f;
    f.end_direction({0,0,0,0,{},PatchOpticalInvalidInput,0,1});
    f.add_hit(one,{0,0});
    f.end_direction({0,0,0,0,{},PatchOpticalQueryError,0,1},PatchQueryUnresolved);
    f.add_hit(one,{0,0});f.vertices.back().basis_x={0,0,1};
    f.end_direction({0,0,0,0,{},PatchOpticalFrameFailure,0,1});
    f.add_hit({{std::numeric_limits<float>::infinity(),0},{0,0}},{0,0});
    f.end_direction({0,0,0,0,{},PatchOpticalInvalidInput,0,1});

    // Source-frame gauge encoded directly in the uploaded input.
    f.add_hit({{.8f,0},{-.6f,0}},{0,0});
    for(auto id:f.patches.back().vertex_indices) f.vertices[id].basis_x={.6f,-.8f,0};
    f.end_direction({1,0,1,0,{1,0,0,0},0,1,0});
    // Known +90-degree rotation y -> z around x, mapping source +z to -y.
    f.add_hit(one,{0,0});
    for(auto id:f.patches.back().vertex_indices)
    { f.vertices[id].direction_drop={0,1,0};f.vertices[id].basis_x={0,0,1}; }
    f.end_direction({1,0,1,0,{1,0,0,0},0,1,0});
    f.add_hit(one,{0,0});
    for(auto id:f.patches.back().vertex_indices) f.vertices[id].direction_drop={0,0,-1};
    f.end_direction({0,0,0,0,{},PatchOpticalFrameFailure,0,1});

    // Random dyadic inputs. Here U,V,fraction have <=8 fractional bits and
    // turns < 2048, so direct weighted q is EXACT even when MSVC long double
    // equals double (<=35 significant bits). No false high-precision claim.
    std::mt19937 rng(0x62bf605u);
    for(unsigned test=0;test<512;++test)
    {
        const float u=float(1u+rng()%254u)/256.0f, v=float(1u+rng()%254u)/256.0f;
        const unsigned count=1u+rng()%7u;
        OpticalExpectation e{};
        std::complex<long double> sum_s{},sum_p{};
        for(unsigned hit=0;hit<count;++hit)
        {
            const float area=float(1u+rng()%16u)/16.0f;
            const float omega=float(1u+rng()%16u)/16.0f;
            f.add_hit({}, {}, area, (test%2u ? -omega : omega),0,
                      test%2u ? PatchCellStatus::RegularNegative : PatchCellStatus::RegularPositive);
            f.hits.back().u=u;f.hits.back().v=v;
            const long double U=u,V=v;
            const long double weights[4]={(1-U)*(1-V),U*(1-V),(1-U)*V,U*V};
            std::complex<long double> s{},p{};long double q=0;
            const auto base=f.patches.back().vertex_indices[0];
            for(unsigned corner=0;corner<4;++corner)
            {
                auto random_field=[&rng](){return float(int(rng()%513u)-256)/128.0f;};
                auto& vertex=f.vertices[base+corner];
                vertex.field={{random_field(),random_field()},{random_field(),random_field()}};
                vertex.optical_cycles={1000+int(rng()%32u),float(rng()%256u)/256.0f};
                s+=weights[corner]*std::complex<long double>{vertex.field.x.real,vertex.field.x.imag};
                p+=weights[corner]*std::complex<long double>{vertex.field.y.real,vertex.field.y.imag};
                q+=weights[corner]*(static_cast<long double>(vertex.optical_cycles.turns)+vertex.optical_cycles.fraction);
                if(test%3u==0u) vertex.optical_cycles.turns+=2000000000;
                else if(test%3u==1u) vertex.optical_cycles.turns-=2000000000;
            }
            const long double amplitude=std::sqrt(static_cast<long double>(area)/omega);
            s*=amplitude;p*=amplitude;
            e.incoherent_s+=static_cast<double>(std::norm(s));e.incoherent_p+=static_cast<double>(std::norm(p));
            constexpr long double two_pi=6.283185307179586476925286766559005768L;
            const long double angle=two_pi*(q-std::floor(q));
            const std::complex<long double> phasor{std::cos(angle),std::sin(angle)};
            sum_s+=s*phasor;sum_p+=p*phasor;
        }
        e.path_s=static_cast<double>(std::norm(sum_s));e.path_p=static_cast<double>(std::norm(sum_p));
        e.field={double(sum_s.real()),double(sum_s.imag()),double(sum_p.real()),double(sum_p.imag())};
        e.evaluated=count;f.end_direction(e,0,test%4u);
    }
    return f;
}

inline void verify_optics_fixture(const OpticsFixture& input, const std::span<const PatchOpticalResult> output)
{
    require_optics(output.size()==input.expected.size(),"Optics output size.");
    for(std::size_t i=0;i<output.size();++i)
    {
        const auto& r=output[i];const auto& e=input.expected[i];
        if(r.flags!=e.flags || r.evaluated_hits!=e.evaluated || r.rejected_hits!=e.rejected
           || r.hit_count!=input.summaries[i].hit_count)
            throw std::runtime_error("Optical flags/counts differ at case "+std::to_string(i));
        if(e.flags & patch_optical_error_mask)
        {
            require_optics(std::isnan(r.regular_partial_incoherent_s)
                && std::isnan(r.regular_partial_path_s)
                && std::isnan(r.regular_partial_path_field.s_real)
                && !r.known_hits_complete(),"Numeric error must not masquerade as zero/valid intensity.");
            continue;
        }
        require_optics(r.known_hits_complete()==(e.flags==0),"Pending/complete distinction.");
        const double power=e.incoherent_s+e.incoherent_p;
        const double field_scale=1.0+std::sqrt(power);
        const auto& a=r.regular_partial_path_field;const auto& b=e.field;
        if(!near_optics(r.regular_partial_incoherent_s,e.incoherent_s)
           || !near_optics(r.regular_partial_incoherent_p,e.incoherent_p)
           || !near_optics(r.regular_partial_path_s,e.path_s,1.0+power)
           || !near_optics(r.regular_partial_path_p,e.path_p,1.0+power)
           || !near_optics(a.s_real,b.s_real,field_scale) || !near_optics(a.s_imag,b.s_imag,field_scale)
           || !near_optics(a.p_real,b.p_real,field_scale) || !near_optics(a.p_imag,b.p_imag,field_scale))
            throw std::runtime_error("Optical values differ at case "+std::to_string(i));
    }
    // Do not let the enormous incoherent sum relax this cancellation test.
    require_optics(output[8].regular_partial_path_field.s_real==1.0
                   && output[8].regular_partial_path_s==1.0,"Compensated cancellation lost the unit wave.");
}

inline void compare_optics_results(const PatchOpticalResult& a,const PatchOpticalResult& b)
{
    require_optics(a.flags==b.flags && a.query_flags==b.query_flags && a.hit_count==b.hit_count
        && a.evaluated_hits==b.evaluated_hits && a.rejected_hits==b.rejected_hits
        && a.first_problem_patch_id==b.first_problem_patch_id
        && a.refinement_hits==b.refinement_hits && a.boundary_hits==b.boundary_hits
        && a.singular_hits==b.singular_hits,"CPU/GPU optical state mismatch.");
    if(a.flags & patch_optical_error_mask) return;
    const double scale=1.0+std::sqrt(b.regular_partial_incoherent_s+b.regular_partial_incoherent_p);
    require_optics(near_optics(a.regular_partial_path_field.s_real,b.regular_partial_path_field.s_real,scale)
        && near_optics(a.regular_partial_path_field.s_imag,b.regular_partial_path_field.s_imag,scale)
        && near_optics(a.regular_partial_path_field.p_real,b.regular_partial_path_field.p_real,scale)
        && near_optics(a.regular_partial_path_field.p_imag,b.regular_partial_path_field.p_imag,scale)
        && near_optics(a.regular_partial_incoherent_s,b.regular_partial_incoherent_s)
        && near_optics(a.regular_partial_incoherent_p,b.regular_partial_incoherent_p)
        && near_optics(a.regular_partial_path_s,b.regular_partial_path_s,scale*scale)
        && near_optics(a.regular_partial_path_p,b.regular_partial_path_p,scale*scale),"CPU/GPU optical values mismatch.");
}
} // namespace rainbow::tests
