#include "wave_optics_test_cases.hpp"
#include <cfenv>
#include <cstdlib>
#include <limits>

namespace
{
using namespace rainbow;using namespace rainbow::tests;
std::size_t checks=0;
void check(bool ok,const char* message){++checks;require_optics(ok,message);}
void phase_tests()
{
    const OpticalField64 f{.25,.5,-1,2};auto g=f;
    for(unsigned i=0;i<4;++i)g=FocalPhase::apply(g,1);
    check(g.s_real==f.s_real&&g.s_imag==f.s_imag&&g.p_real==f.p_real&&g.p_imag==f.p_imag,"Four quarters not identity.");
    auto fixture=make_wave_fixture();auto p=fixture.optical.view();
    for(unsigned j=0;j<2;++j)
    {
        FocalPhaseEstimate e{};
        check(FocalPhase::estimate(p.vertices,p.vertex_count,fixture.optical.patches[j],p.incident_direction,fixture.focal,e)==0,"Focal estimate failed.");
        check(e.extra_quarter_turn==j,"Focal derivative sign incorrect.");
        check(wave_near(e.radial_derivative,j?.1:-.1,1e-6),"Bilinear radial derivative reference.");
    }
    for(unsigned q=0;q<4;++q)
    {
        auto config=fixture.focal;config.quarter_turn_offsets[1]=q;
        FocalOpticalResult out{};
        const auto raw=PatchOpticalEvaluator::evaluate_direction(p,0,&config,&out);
        const double expected[4]={2,4,2,0};
        check(out.valid()&&out.corrected_hits==2,"Focal evaluator validity.");
        check(wave_near(out.intensity_s,expected[q],2e-6),"Two-wave quarter-phase interference.");
        const auto plain=PatchOpticalEvaluator::evaluate_direction(p,0);
        check(raw.regular_partial_path_s==plain.regular_partial_path_s
            &&raw.regular_partial_incoherent_s==plain.regular_partial_incoherent_s,"Original path output changed.");
    }
    FocalOpticalResult a{},b{};
    static_cast<void>(PatchOpticalEvaluator::evaluate_direction(p,0,&fixture.focal,&a)); // checked via outputs
    auto shifted=fixture.focal;for(auto& offset:shifted.quarter_turn_offsets)offset=3;
    static_cast<void>(PatchOpticalEvaluator::evaluate_direction(p,0,&shifted,&b));
    check(a.intensity_s==b.intensity_s,"Common phase offset changed intensity.");
    fixture.optical.hits[0].flags=BilinearBoundary;
    auto raw=PatchOpticalEvaluator::evaluate_direction(p,0,&fixture.focal,&b);
    check((raw.flags&PatchOpticalBoundaryPending) && (b.flags&FocalInputPending)
          &&std::isnan(b.intensity_s),"Pending geometry hidden by focal stage.");
    auto bad=fixture.focal;bad.grid_width=0;
    raw=PatchOpticalEvaluator::evaluate_direction(p,0,&bad,&b);
    check((b.flags&FocalInvalidGeometry)!=0,"Invalid layout accepted.");
    auto patch=fixture.optical.patches[1];patch.vertex_indices[0]++;
    FocalPhaseEstimate e{};e.radial_derivative=123;
    check(FocalPhase::estimate(p.vertices,p.vertex_count,patch,p.incident_direction,fixture.focal,e)!=0
          &&e.radial_derivative==123,"Failure modified focal estimate.");
    // Constant outgoing angle has undefined sign, never assign an arbitrary quarter.
    fixture.optical.hits[0].flags=0;
    for(auto& v:fixture.optical.vertices)v.direction_drop=fixture.optical.directions[0];
    check(FocalPhase::estimate(p.vertices,p.vertex_count,fixture.optical.patches[1],p.incident_direction,fixture.focal,e)
          ==FocalDerivativePending,"Zero derivative was silently resolved.");
    // An emitting cell centered exactly on the axis: no radial direction.
    FocalPhaseConfig centered{2,2,1,0,{0,0,0,0}};
    OutgoingPatch origin{};origin.vertex_indices[0]=0;origin.vertex_indices[1]=1;origin.vertex_indices[2]=2;origin.vertex_indices[3]=3;
    check(FocalPhase::estimate(p.vertices,p.vertex_count,origin,p.incident_direction,centered,e)
          ==FocalDerivativePending,"Axis-centered radial derivative accepted.");
}
void diffraction_tests()
{
    double sigma=999;
    constexpr float r[10]={.1f,.2f,.3f,.4f,.5f,.6f,.7f,.8f,.9f,1.f};
    constexpr double s[10]={.70,.45,.30,.25,.22,.20,.18,.17,.16,.15};
    for(unsigned i=0;i<10;++i)check(RainbowDiffraction::table_sigma_degrees(r[i],sigma)&&sigma==s[i],"Table II mismatch.");
    sigma=999;check(!RainbowDiffraction::table_sigma_degrees(1.01f,sigma)&&sigma==999,"Table extrapolation.");
    check(!RainbowDiffraction::table_sigma_degrees(.01f,sigma),"Table lower extrapolation.");
    check(RainbowDiffraction::table_sigma_degrees(.35f,sigma)&&wave_near(sigma,.275,1e-7),"Radius interpolation.");
    auto g=make_diffraction_fixture();auto p=g.view();
    unsigned primary=0,secondary=0;
    for(auto e:g.edges){primary+=(e.flags&TransitionPrimary)!=0;secondary+=(e.flags&TransitionSecondary)!=0;}
    check(primary==g.columns&&secondary==g.columns,"Rainbow labels missing.");
    bool saw_primary=false,saw_secondary=false,saw_taper=false;
    for(std::uint32_t i=0;i<g.input.size();i+=7)
    {
        const auto out=RainbowDiffraction::filter(p,i);
        check(out.valid(),"Valid fixture filtered to invalid.");
        if(out.flags&DiffractionFiltered)
        {
            saw_primary|=out.transition_kind==TransitionPrimary;saw_secondary|=out.transition_kind==TransitionSecondary;
            saw_taper|=out.blend>0&&out.blend<1;
            check(out.sigma_rad==p.config.primary_sigma_rad*(out.transition_kind==TransitionSecondary?2:1),"Secondary width is not doubled.");
        }
    }
    check(saw_primary&&saw_secondary&&saw_taper,"Missing primary/secondary/taper test coverage.");
    // Independent all-grid convolution: does NOT reuse the cap-window search.
    for(std::uint32_t row:{126u,129u,137u,140u})for(std::uint32_t col:{0u,1u,71u})
    {
        const auto id=row*g.columns+col;const auto out=RainbowDiffraction::filter(p,id);
        check((out.flags&DiffractionFiltered)!=0,"Reference target not filtered.");
        const long double dt=RainbowDiffraction::pi/g.rows,dp=2*RainbowDiffraction::pi/g.columns;
        const long double theta=(row+.5L)*dt;long double den=0,num_s=0,num_p=0;
        for(std::uint32_t rr=0;rr<g.rows;++rr)for(std::uint32_t cc=0;cc<g.columns;++cc)
        {
            const long double t=(rr+.5L)*dt,dphi=(static_cast<long double>(cc)-col)*dp;
            const long double dot=std::cos(theta)*std::cos(t)+std::sin(theta)*std::sin(t)*std::cos(dphi);
            const long double d=std::acos(std::clamp(dot,-1.0L,1.0L));
            if(d>p.config.support_sigma*out.sigma_rad+RainbowDiffraction::angular_boundary_tolerance)continue;
            const long double w=std::exp(-.5L*(d/out.sigma_rad)*(d/out.sigma_rad))
                *2*std::sin(t)*std::sin(dt/2)*dp;
            const auto& x=g.input[rr*g.columns+cc];den+=w;num_s+=w*x.intensity_s;num_p+=w*x.intensity_p;
        }
        check(wave_near(out.intensity_s,double(g.input[id].intensity_s+out.blend*(num_s/den-g.input[id].intensity_s)),2e-10),"Spherical stencil differs from full sum.");
        check(wave_near(out.intensity_p,double(g.input[id].intensity_p+out.blend*(num_p/den-g.input[id].intensity_p)),2e-10),"Spherical p stencil differs from full sum.");
    }
    // Constant preservation with independently supplied transition labels.
    for(auto& x:g.input){x.intensity_s=3.25;x.intensity_p=.125;}
    auto value=RainbowDiffraction::filter(p,137*g.columns);
    check(wave_near(value.intensity_s,3.25,1e-14)&&wave_near(value.intensity_p,.125,1e-14),"Gaussian did not preserve constant.");
    // No hole filling, even if the omitted field could be renormalized away.
    g.input[137*g.columns+1].flags=FocalInputPending;
    value=RainbowDiffraction::filter(p,137*g.columns);
    check(!value.valid()&&(value.flags&DiffractionStencilUnavailable)&&std::isnan(value.intensity_s),"Invalid stencil silently filled.");
    g.input[137*g.columns+1].flags=0;
    g.edges[137*g.columns+1].flags=TransitionUnknown;
    value=RainbowDiffraction::filter(p,137*g.columns);
    check(!value.valid()&&(value.flags&DiffractionDetectionUnavailable),"Unknown transition silently ignored.");
    for(auto& e:g.edges)e.flags=0;
    value=RainbowDiffraction::filter(p,10*g.columns);
    check(value.flags==0&&value.intensity_s==3.25,"No-edge region changed.");
    for(auto& x:g.input){x.intensity_s=0;x.intensity_p=0;}
    value=RainbowDiffraction::filter(p,10*g.columns);
    check(value.valid()&&value.intensity_s==0,"Physical zero replaced by NaN.");
    p.config.primary_sigma_rad=std::numeric_limits<double>::quiet_NaN();
    check(!RainbowDiffraction::valid_config(p),"NaN config accepted.");
    // One-pixel kernel at coarse resolution: warn, do not widen sigma.
    g=make_diffraction_fixture();p=g.view();p.config.primary_sigma_rad=.15*RainbowDiffraction::pi/180;
    value=RainbowDiffraction::filter(p,137*g.columns);
    check(value.valid()&&(value.flags&DiffractionUnderresolved),"Coarse kernel unflagged.");
}
}
int main()
{
    try
    {
        require_optics(std::fegetround()==FE_TONEAREST,"Round-to-nearest required.");
        phase_tests();diffraction_tests();
        std::cout<<"Wave optics CPU: "<<checks<<" checks passed. long-double digits="
            <<std::numeric_limits<long double>::digits<<" (MSVC may equal double).\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return EXIT_FAILURE;}
}
