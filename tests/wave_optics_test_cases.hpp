#pragma once
#include "patch_optics_test_cases.hpp"
#include <rainbow/focal_phase.hpp>
#include <rainbow/rainbow_diffraction.hpp>
#include <algorithm>
#include <iostream>

namespace rainbow::tests
{
inline bool wave_near(double a,double b,double tolerance=2e-6)
{return std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=tolerance*(1+std::abs(b));}
inline void compare_focal(const FocalOpticalResult& a,const FocalOpticalResult& b)
{
    require_optics(a.flags==b.flags && a.corrected_hits==b.corrected_hits
        && a.extra_quarter_turn_hits==b.extra_quarter_turn_hits
        && a.first_problem_patch_id==b.first_problem_patch_id,"Focal integer result mismatch.");
    for(unsigned f=0;f<4;++f)require_optics(a.family_hits[f]==b.family_hits[f]
        && wave_near(a.family_incoherent[f],b.family_incoherent[f]),"Focal family mismatch.");
    if(a.valid())
    {
        require_optics(wave_near(a.intensity_s,b.intensity_s)&&wave_near(a.intensity_p,b.intensity_p),"Focal intensity mismatch.");
        require_optics(wave_near(a.field.s_real,b.field.s_real)&&wave_near(a.field.s_imag,b.field.s_imag)
            &&wave_near(a.field.p_real,b.field.p_real)&&wave_near(a.field.p_imag,b.field.p_imag),"Focal field mismatch.");
    }
    else require_optics(std::isnan(a.intensity_s)&&std::isnan(a.intensity_p),"Unavailable focal intensity must be NaN.");
}
inline void compare_diffraction(const DiffractionResult& a,const DiffractionResult& b)
{
    require_optics(a.flags==b.flags && a.transition_kind==b.transition_kind,"Diffraction flags mismatch.");
    require_optics(wave_near(a.sigma_rad,b.sigma_rad,1e-11)&&wave_near(a.blend,b.blend,1e-9),"Diffraction selection mismatch.");
    if(a.valid())require_optics(wave_near(a.intensity_s,b.intensity_s,1e-9)
        &&wave_near(a.intensity_p,b.intensity_p,1e-9),"Diffraction intensity mismatch.");
    else require_optics(std::isnan(a.intensity_s)&&std::isnan(a.intensity_p),"Unavailable diffraction must be NaN.");
}
struct WaveFixture
{
    OpticsFixture optical;
    FocalPhaseConfig focal{5,5,1.0f,0,{0,0,0,0}};
};
inline WaveFixture make_wave_fixture(std::uint32_t count=129)
{
    WaveFixture x;auto& f=x.optical;
    f.incident={0,0,1};f.incident_basis={1,0,0};
    f.vertices.resize(100);
    for(unsigned family=0;family<4;++family)
        for(unsigned iy=0;iy<5;++iy)for(unsigned ix=0;ix<5;++ix)
        {
            const double U=-1+0.5*ix,V=-1+0.5*iy;
            const double beta=family==1?-0.1:0.1;
            const double t=0.65+beta*(U*U+V*V);
            auto& v=f.vertices[family*25+iy*5+ix];
            v.status=VertexStatus::Valid;v.direction_drop=Vec3{float(std::sin(t)),0,float(std::cos(t))}.normalized();
            v.basis_x={0,1,0};v.field={{1,0},{0,0}};v.optical_cycles={1000,0};
        }
    // Arithmetic fixture only: assigned areas and hit coords, not a new GAS model.
    for(unsigned family:{1u,2u})
    {
        OutgoingPatch p{};p.patch_id=family*16+11;p.status=PatchCellStatus::RegularPositive;
        p.signed_solid_angle_sr=1;p.incident_area_drop2=1;
        const auto v=family*25+13;
        p.vertex_indices[0]=v;p.vertex_indices[1]=v+1;p.vertex_indices[2]=v+5;p.vertex_indices[3]=v+6;
        f.patches.push_back(p);
    }
    for(std::uint32_t i=0;i<count;++i)
    {
        f.directions.push_back(Vec3{float(std::sin(.65)),0,float(std::cos(.65))}.normalized());
        f.hits.push_back({0,27,0,0,.5f,.5f,1,0});f.hits.push_back({1,43,0,0,.5f,.5f,1,0});
        PatchQuerySummary s{};s.hit_count=2;s.candidate_count=2;s.regular_hits=2;f.summaries.push_back(s);
        f.offsets.push_back(f.hits.size());
    }
    return x;
}
struct DiffractionFixture
{
    std::uint32_t rows=180,columns=72;
    DiffractionConfig config{2*RainbowDiffraction::pi/180,4,0.25,3};
    std::vector<FocalOpticalResult> input;
    std::vector<RainbowTransition> edges;
    [[nodiscard]] DiffractionParams view() noexcept
    {return {input.data(),edges.data(),nullptr,rows,columns,config};}
};
inline DiffractionFixture make_diffraction_fixture(std::uint32_t rows=180,std::uint32_t columns=72)
{
    DiffractionFixture g;g.rows=rows;g.columns=columns;
    g.input.resize(std::size_t(rows)*columns);g.edges.resize(g.input.size());
    for(std::uint32_t r=0;r<rows;++r)for(std::uint32_t c=0;c<columns;++c)
    {
        const double degrees=180*(double(r)+.5)/rows;
        const double phi=2*RainbowDiffraction::pi*c/columns;
        const bool primary=degrees>138+0.8*std::sin(phi),secondary=degrees<129+0.5*std::sin(phi);
        const double I=1+(primary?8:0)+(secondary?2:0);
        auto& v=g.input[std::size_t(r)*columns+c];
        v.intensity_s=.8*I;v.intensity_p=.2*I;
        v.family_hits[2]=primary?2u:0u;v.family_hits[3]=secondary?2u:0u;
        v.family_incoherent[2]=primary?8:0;v.family_incoherent[3]=secondary?2:0;
    }
    auto p=g.view();for(std::uint32_t i=0;i<g.input.size();++i)g.edges[i]=RainbowDiffraction::detect(p,i);
    return g;
}
} // namespace rainbow::tests
