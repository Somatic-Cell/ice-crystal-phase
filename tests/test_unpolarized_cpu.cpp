#include "unpolarized_test_support.hpp"
#include <rainbow/rainbow_diffraction.hpp>
#include <array>
#include <rainbow/patch_failure_report.hpp>
#include <chrono>
#include <fstream>
#include <iterator>

using namespace unpolarized_tests;
namespace
{
void interfaces()
{
    auto j=JonesResponse32::identity();
    require(near(j.apply_to({{.3f,.2f},{-.4f,.7f}}),Field32{{.3f,.2f},{-.4f,.7f}}),"Jones identity failed.");
    unsigned tir=0,transmitted=0;
    for(unsigned interior=0;interior<2;++interior)for(unsigned angle=0;angle<90;++angle)
    {
        const float a=float(angle*QueryDirectionGrid::pi/180.0);
        const Vec3 w=Vec3{std::sin(a),0,-std::cos(a)}.normalized();
        const auto f=TransverseFrame::from_direction(w);
        j.column[0]={{.3f,.2f},{-.7f,.1f}};j.column[1]={{.5f,-.8f},{.9f,.4f}};
        const float ni=interior?1.3314f:1.0f,nt=interior?1.0f:1.3314f;
        const auto both=DielectricInterface::split(JonesResponseRay{w,f,j},{0,0,1},ni,nt);
        require(both.is_valid,"Test interface not valid.");
        tir+=!both.has_transmission;transmitted+=both.has_transmission;
        for(unsigned k=0;k<2;++k)
        {
            const auto one=DielectricInterface::split(PolarizedRay{w,f,j.column[k]},{0,0,1},ni,nt);
            require(one.is_valid==both.is_valid&&one.has_transmission==both.has_transmission&&one.diagnostics==both.diagnostics,"Shared interface state differs.");
            require(std::memcmp(&one.reflected.direction,&both.reflected.direction,sizeof(Vec3))==0,"Reflected directions differ.");
            require(near(one.reflected.field,both.reflected.field.column[k]),"Reflected matrix column differs.");
            if(one.has_transmission)
            {
                require(std::memcmp(&one.transmitted.direction,&both.transmitted.direction,sizeof(Vec3))==0,"Transmitted directions differ.");
                require(near(one.transmitted.field,both.transmitted.field.column[k]),"Transmitted matrix column differs.");
            }
        }
    }
    require(tir>0&&transmitted>0,"TIR/transmission not covered.");
}
void algebra()
{
    Traced a;a.config=settings(.4f,2,true).make_config();a.vertices.resize(8);a.second.resize(8);
    Geometry g;g.patches.resize(2);g.lookup.assign(2,0xffffffffu);
    const Vec3 direction=Vec3{-.4f,.5f,.7f}.normalized();
    TransverseFrame frame{};require(PolarizationTransport::try_make_scattering_frame(a.config.incident_direction,a.config.incident_basis_x,direction,frame),"Bad target frame.");
    Query q;q.direction=direction;
    for(unsigned p=0;p<2;++p)
    {
        auto& patch=g.patches[p];patch.status=PatchCellStatus::RegularPositive;
        patch.patch_id=p;patch.incident_area_drop2=patch.signed_solid_angle_sr=1;
        for(unsigned k=0;k<4;++k)
        {
            unsigned i=4*p+k;patch.vertex_indices[k]=i;
            auto& v=a.vertices[i];v.status=VertexStatus::Valid;v.direction_drop=direction;v.basis_x=frame.e0;
            v.field={{1,0},{0,0}};a.second[i]={{0,0},{1,0}};
            v.optical_cycles={0,p?.5f:0.0f};
        }
        q.hits.push_back({p,p,0,0,.3f,.4f,1,0});
    }
    q.summary.hit_count=q.summary.candidate_count=2;q.summary.regular_hits=2;q.offsets[1]=2;
    auto params=q.params(a,g,true);
    auto result=PatchOpticalEvaluator::evaluate_direction(params,0);
    require(result.known_hits_complete(),"Algebra input failed.");
    require(near(result.regular_partial_incoherent_s+result.regular_partial_incoherent_p,2,1e-6),"Unit intensity/state weights incorrect.");
    require(near(result.regular_partial_path_s+result.regular_partial_path_p,0,1e-6),"Same-column interference lost.");
    // One rank-one response: both input columns map onto the same output axis.
    // Columns MUST NOT interfere with each other.
    q.hits.resize(1);q.summary.hit_count=q.summary.candidate_count=q.summary.regular_hits=1;q.offsets[1]=1;
    for(unsigned k=0;k<4;++k)a.second[k]={{1,0},{0,0}};
    params=q.params(a,g,true);result=PatchOpticalEvaluator::evaluate_direction(params,0);
    require(near(result.regular_partial_path_s+result.regular_partial_path_p,1,1e-6),"Independent input columns interfered.");
    for(unsigned k=0;k<4;++k)a.second[k]={};
    result=PatchOpticalEvaluator::evaluate_direction(params,0);
    require(near(result.regular_partial_path_s+result.regular_partial_path_p,.5,1e-6),"Ideal polarizer did not transmit half.");
    a.second[2].x.real=std::numeric_limits<float>::quiet_NaN();
    result=PatchOpticalEvaluator::evaluate_direction(params,0);
    require((result.flags&patch_optical_error_mask)!=0&&!std::isfinite(result.regular_partial_path_s),"Bad second column was hidden.");
    params.second_input_fields=nullptr;
    require(PatchOpticalEvaluator::evaluate_direction(params,0).flags&PatchOpticalInvalidInput,"Missing paired buffer silently fell back.");
    params=q.params(a,g,true);--params.second_input_count;
    require(PatchOpticalEvaluator::evaluate_direction(params,0).flags&PatchOpticalInvalidInput,"Mismatched paired count accepted.");
}
void serialization(const std::filesystem::path& requested)
{
    PatchFailureSnapshot snap{};
    snap.input_polarization=IncidentPolarization::Unpolarized;
    snap.config=settings(.4f,2,true).make_config();snap.theta_count=4;snap.phi_count=8;snap.stored_patch_count=1;
    FailedPatchWitness w{};w.patch.status=PatchCellStatus::RegularPositive;
    w.patch.incident_area_drop2=1;
    const Vec3 directions[4]={Vec3{-.1f,-.1f,1}.normalized(),Vec3{.1f,-.1f,1}.normalized(),
        Vec3{-.1f,.1f,1}.normalized(),Vec3{.1f,.1f,1}.normalized()};
    BilinearPatchGeometry g{};
    for(unsigned k=0;k<4;++k)
    {
        auto& v=w.vertices[k];v.status=VertexStatus::Valid;v.direction_drop=directions[k];
        v.basis_x=TransverseFrame::from_direction(directions[k]).e0;v.field={{1,0},{0,0}};
        w.second_input_fields[k]={{0,0},{1,0}};w.patch.vertex_indices[k]=k;g.corners[k]=directions[k];
    }
    w.patch.signed_solid_angle_sr=g.signed_solid_angle();snap.patches.push_back(w);
    FailedDirectionWitness d{};d.direction_id=0;d.direction={0,0,1};
    d.query.hit_count=d.query.candidate_count=d.query.regular_hits=1;
    d.optical.hit_count=d.optical.evaluated_hits=1;
    d.hits.push_back({0,0,0,0,.5f,.5f,1,0});
    d.focal.flags=FocalDerivativePending;d.focal.intensity_s=d.focal.intensity_p=patch_optical_nan;
    d.diffraction.flags=DiffractionInputUnavailable;
    snap.directions.push_back(d);snap.focal_unavailable_directions=1;snap.diffraction_unavailable_directions=1;snap.nonfinite_directions=1;
    const auto path=requested.empty()?std::filesystem::temp_directory_path()/
        ("rainbow_unpolarized_cpu_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".json"):requested;
    PatchFailureReport(std::move(snap)).write_json(path);
    std::string data;
    {std::ifstream file(path);data.assign(std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>());}
    require(data.find("rainbow_patch_failure_witness_v2")!=std::string::npos,"Paired report version lost.");
    require(data.find("second_input_field_bits")!=std::string::npos,"Paired report lacks column 1.");
    require(data.find("incident_field_bits")==std::string::npos,"Paired report mislabeled a single input.");
    require(data.find("[0.5,0,0,0.5]")!=std::string::npos,"Paired report input statistics incorrect.");
    if(requested.empty())std::filesystem::remove(path);
}
void optical_case(const RaindropSettings& s)
{
    auto pair=trace(s,true);auto sx=s,sy=s;sx.incident_field={{1,0},{0,0}};sy.incident_field={{0,0},{1,0}};
    const auto x=trace(sx,false),y=trace(sy,false);check_trace(pair,x,y);
    auto g=build(pair);
    const auto directions=QueryDirectionGrid::make(pair.config,90,180).directions();
    // The 36 saved diagnostic directions, including the 12 formerly pending
    // directions, plus ordinary forward/side/backward directions. The new
    // responses below are retraced, not manufactured from the old x snapshot.
    const std::vector<unsigned> selected={73,1758,4220,8004,11172,11173,11236,11237,11352,11353,11362,11407,11416,11417,11532,11533,11542,11546,11583,11587,11596,11597,11628,11630,11679,11681,11722,11726,11763,11767,11808,11810,11859,11861,11906,11943,11988,11990,12039,12041,14999};
    unsigned complete=0,folded_hits=0;
    auto rotated=pair;
    // Postmultiply every response by one unitary matrix. Geometry and the grid
    // stay exactly unchanged. Complex circular columns exercise cross-coupling.
    constexpr float h=.7071067811865475244f;
    for(unsigned i=0;i<rotated.vertices.size();++i)
    {
        JonesResponse32 j{{pair.vertices[i].field,pair.second[i]}};
        rotated.vertices[i].field=j.apply_to({{h,0},{0,h}});
        rotated.second[i]=j.apply_to({{h,0},{0,-h}});
    }
    for(auto i:selected)
    {
        auto q=query(pair,g,directions.at(i));
        FocalOpticalResult f{},fx{},fy{},fr{};
        const auto p=PatchOpticalEvaluator::evaluate_direction(q.params(pair,g,true),0,&g.focal,&f,g.view());
        const auto px=PatchOpticalEvaluator::evaluate_direction(q.params(x,g,false),0,&g.focal,&fx,g.view());
        const auto py=PatchOpticalEvaluator::evaluate_direction(q.params(y,g,false),0,&g.focal,&fy,g.view());
        check_average(p,f,px,fx,py,fy);
        const auto pr=PatchOpticalEvaluator::evaluate_direction(q.params(rotated,g,true),0,&g.focal,&fr,g.view());
        require(pr.flags==p.flags&&fr.flags==f.flags,"Unitary input change changed availability.");
        if(!p.known_hits_complete()||!f.valid())continue;
        ++complete;folded_hits+=p.folded_evaluated_hits();
        require(near(pr.regular_partial_incoherent_s,p.regular_partial_incoherent_s)&&near(pr.regular_partial_incoherent_p,p.regular_partial_incoherent_p),"Unpolarized incoherent not input-basis invariant.");
        require(near(pr.regular_partial_path_s,p.regular_partial_path_s,2e-5)&&near(pr.regular_partial_path_p,p.regular_partial_path_p,2e-5),"Unpolarized path not input-basis invariant.");
        require(near(fr.intensity_s,f.intensity_s,2e-5)&&near(fr.intensity_p,f.intensity_p,2e-5),"Unpolarized focal not input-basis invariant.");
    }
    require(complete>selected.size()/2,"Too few valid optical comparisons.");
    if(!s.force_sphere&&s.grid_width==129)require(folded_hits==12,"Nonspherical folded response not exercised.");
    std::cout<<"pair radius="<<s.radius_mm<<" grid="<<s.grid_width<<": "<<pair.calls<<" geometry queries (one input pass); "
             <<complete<<" complete directions; "<<folded_hits<<" folded contributions; two-state mean and unitary invariance passed.\n";
}
}
int main(int argc,char** argv)
{
    try
    {
        if(argc>2)throw std::invalid_argument("Optional argument: new synthetic report path");
        serialization(argc==2?std::filesystem::path(argv[1]):std::filesystem::path{});
        interfaces();algebra();optical_case(settings(.4f,17,true));optical_case(settings(1.0f,129));
        std::cout<<checks<<" unpolarized checks passed (CPU).\n";return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
