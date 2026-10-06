#include <rainbow/patch_failure_report.hpp>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
using namespace rainbow;
unsigned checks=0;
void require(const bool condition,const char* message)
{++checks;if(!condition)throw std::runtime_error(message);}
BilinearPatchGeometry regular()
{return {{{-.2f,-.1f,1.f},{.2f,-.1f,1.f},{-.2f,.1f,1.f},{.2f,.1f,1.f}}};}
BilinearPatchGeometry thin()
{return {{{0,0,1},{0x1p-10f,0x1p-10f,1},{0x1p-10f,0x1p-10f+0x1p-32f,1},{0x1p-9f,0x1p-9f+0x1p-32f,1}}};}
BilinearPatchGeometry fold()
{auto g=regular();std::swap(g.corners[2],g.corners[3]);return g;}
FailedPatchWitness witness(const BilinearPatchGeometry& g,const std::uint32_t index,const RaindropTraceConfig& config)
{
    FailedPatchWitness w{};w.compact_index=index;w.patch.patch_id=2*index;
    PatchBuildLayout layout{};
    require(PatchBuildLayout::try_make(config,layout),"Test layout failed");
    layout.corner_indices(w.patch.patch_id,w.patch.vertex_indices);
    w.patch.incident_area_drop2=layout.incident_area_drop2;
    const int orientation=g.certified_orientation();
    w.patch.status=orientation>0?PatchCellStatus::RegularPositive:
                   orientation<0?PatchCellStatus::RegularNegative:PatchCellStatus::NeedsRefinement;
    w.patch.signed_solid_angle_sr=orientation?g.signed_solid_angle():std::numeric_limits<float>::quiet_NaN();
    for(unsigned k=0;k<4;++k)
    {
        w.vertices[k].direction_drop=g.corners[k];w.vertices[k].status=VertexStatus::Valid;
        w.vertices[k].position_drop={0,0,0};w.vertices[k].basis_x={1,0,0};
        w.vertices[k].field={{1,0},{0,0}};w.vertices[k].optical_cycles={123,0.25f};
    }
    return w;
}
void numerical_tests()
{
    auto g=regular();const auto a=PatchRegularityAudit::inspect(g);
    require(a.finite_input,"Regular input marked invalid");
    require(a.host_orientation32==1 && a.interval_orientation64==1,"Regular sign");
    require(a.signed_omega64>0,"Regular omega");
    std::swap(g.corners[1],g.corners[2]);const auto n=PatchRegularityAudit::inspect(g);
    require(n.host_orientation32==-1 && n.interval_orientation64==-1,"Negative sign");
    const auto thin_audit=PatchRegularityAudit::inspect(thin());
    require(thin_audit.host_orientation32==0 && thin_audit.interval_orientation64==1,"FP64 must certify this thin regular fixture");
    g=fold();const auto f=PatchRegularityAudit::inspect(g);
    require(f.host_orientation32==0 && f.interval_orientation64==0,"Fold must stay uncertified");
    require(f.jacobian_value64[0]>0 && f.jacobian_value64[2]<0,"Fold sign witness");
    g=regular();g.corners[1]=g.corners[0];g.corners[3]=g.corners[2];
    const auto z=PatchRegularityAudit::inspect(g);
    require(z.interval_orientation64==0,"Exact degeneracy must stay unresolved");
    g=regular();g.corners[0].x=std::numeric_limits<float>::quiet_NaN();
    require(!PatchRegularityAudit::inspect(g).finite_input,"NaN input must be recorded");
    // Independent direct evaluation of J inside bilinear parameter domain vs
    // its four-corner bilinear interpolation. This is a floating-point check;
    // Python tests separately verify exact signs with rational arithmetic.
    std::uint32_t random=177;
    for(unsigned i=0;i<256;++i)
    {
        g=regular();
        for(auto& p:g.corners)
        {
            random=1664525u*random+1013904223u;p.x+=(float(random&1023u)-511.f)*0.0002f;
            random=1664525u*random+1013904223u;p.y+=(float(random&1023u)-511.f)*0.0002f;
            random=1664525u*random+1013904223u;p.z+=(float(random&1023u)-511.f)*0.0001f;
        }
        const auto r=PatchRegularityAudit::inspect(g);
        const double u=.37,v=.61;
        const auto& j=r.jacobian_value64;
        const double b=(1-u)*(1-v)*j[0]+u*(1-v)*j[1]+(1-u)*v*j[2]+u*v*j[3];
        require(std::abs(b-PatchRegularityAudit::jacobian_at(g,u,v))<2e-15,"J must be bilinear");
        for(unsigned k=0;k<4;++k)
            require(r.jacobian64[k].lo<=j[k] && j[k]<=r.jacobian64[k].hi,"Double central value enclosure");
    }
}
PatchFailureSnapshot fixtures()
{
    PatchFailureSnapshot s{};
    s.config.grid_width=4097;s.config.grid_height=2;s.config.grid_half_extent=1.01f;
    s.theta_count=2;s.phi_count=2;s.config.radius_mm=1;s.stored_patch_count=132;
    std::uint32_t state=1541;
    for(std::uint32_t i=0;i<132;++i)
    {
        auto g=regular();
        if(i==1)g=fold();
        else if(i==2){g.corners[1]=g.corners[0];g.corners[3]=g.corners[2];}
        else if(i==3)std::swap(g.corners[1],g.corners[2]);
        else if(i==4)g=thin();
        else if(i>=5)
        {
            for(auto& p:g.corners)
            {
                state=1664525u*state+1013904223u;p.x+=(float(state&1023u)-511.f)*0.0005f;
                state=1664525u*state+1013904223u;p.y+=(float(state&1023u)-511.f)*0.0005f;
                p=p.normalized();
            }
        }
        s.patches.push_back(witness(g,i,s.config));
    }
    FailedDirectionWitness d{};d.direction_id=0;d.direction={0,0,1};
    d.query.hit_count=d.query.candidate_count=1;d.query.refinement_hits=1;d.query.flags=PatchQueryRefinementHit;
    d.optical.hit_count=d.optical.rejected_hits=d.optical.refinement_hits=1;
    d.optical.flags=PatchOpticalRefinementPending;d.optical.first_problem_patch_id=2;
    d.focal.flags=FocalInputPending;d.focal.first_problem_patch_id=2;
    d.focal.field={patch_optical_nan,patch_optical_nan,patch_optical_nan,patch_optical_nan};
    d.focal.intensity_s=d.focal.intensity_p=patch_optical_nan;
    d.diffraction.flags=DiffractionInputUnavailable;
    d.diffraction.intensity_s=d.diffraction.intensity_p=patch_optical_nan;
    d.hits.push_back({1,2,0,0,.25f,.75f,1.f,0.f});
    s.directions.push_back(d);s.optical_incomplete_directions=1;s.focal_unavailable_directions=1;
    s.diffraction_unavailable_directions=1;s.nonfinite_directions=1;
    return s;
}
}
int main(int argc,char** argv)
{
    try
    {
        numerical_tests();
        const bool keep=argc==2;
        const auto path=keep?std::filesystem::path(argv[1]):std::filesystem::temp_directory_path()/
            ("rainbow_patch_report_test_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".json");
        const PatchFailureReport report(fixtures());report.write_json(path);
        
        std::ifstream in(path,std::ios::binary);
        const std::string text((std::istreambuf_iterator<char>(in)),{});
        in.close();

        require(text.find("\"report_complete\":true")!=std::string::npos,"Missing completion marker");
        require(text.find("\"origin\":\"synthetic_test\"")!=std::string::npos,"Fixture provenance lost");
        require(text.find("\"focal_intensity\":[null,null]")!=std::string::npos,"NaN must be valid JSON null");
        bool refused=false;try{report.write_json(path);}catch(const std::invalid_argument&){refused=true;}
        require(refused,"Existing report must not be overwritten");
        auto broken=fixtures();broken.directions[0].hits[0].patch_id=17;
        const PatchFailureReport bad(std::move(broken));
        auto impossible=path;impossible+=".invalid";
        refused=false;try{bad.write_json(impossible);}catch(const std::invalid_argument&){refused=true;}
        require(refused && !std::filesystem::exists(impossible),"Inconsistent witness must be rejected before writing");
        if(!keep)std::filesystem::remove(path);
        std::cout<<checks<<" checks passed (CPU diagnostic only).\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
