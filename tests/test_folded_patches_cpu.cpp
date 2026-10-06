#include "folded_test_support.hpp"
#include <rainbow/rainbow_diffraction.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace
{
using namespace rainbow;using namespace folded_tests;
void component_tests()
{
    FoldedPatchGeometry g{};g.j[0]=1;g.j[1]=-1;g.j[2]=-1;g.j[3]=1;g.scale=1;
    FoldedPatchRecord r{};
    require(g.components(r)&&r.branch_count==4,"Disconnected same-sign components were merged.");
    double area=0;
    for(unsigned i=0;i<r.branch_count;++i)
    {
        const auto& b=r.branches[i];double lo=0,hi=0;
        require(g.section(g.j,b.sign,(b.u0+b.u1)*.5,lo,hi),"Empty component.");
        area+=(b.u1-b.u0)*(hi-lo);
    }
    require(near(area,1),"Checkerboard domains do not partition the square.");
}
void double_root_test()
{
    // A symmetric genuine folded chord patch; two distinct UV roots, one on
    // each side of the diagonal. No triangle substitution or nearest-hit rule.
    Vec3 corners[4]={{0,0,1},Vec3{.2f,0,1}.normalized(),Vec3{.2f,0,1}.normalized(),Vec3{.4f,.1f,1}.normalized()};
    const auto record=FoldedPatchBuilder::geometry(corners,0,0);
    require(record.flags==0&&record.branch_count==2,"Synthetic fold preparation failed.");
    require(near(record.branches[0].area_fraction,.5,1e-8)&&near(record.branches[1].area_fraction,.5,1e-8),"Symmetric branch area is not one half.");
    BilinearPatchGeometry geometry{};for(unsigned i=0;i<4;++i)geometry.corners[i]=corners[i];
    const Vec3 direction=geometry.evaluate(.23f,.71f).normalized();
    const auto hits=BilinearPatchIntersector::intersect(geometry,direction);
    require(hits.count==2 && !(hits.flags&(BilinearUnresolved|BilinearInvalidInput)),"Two folded roots were not preserved.");
    const int a=FoldedPatchGeometry::branch_at(record,hits.hits[0].u,hits.hits[0].v);
    const int b=FoldedPatchGeometry::branch_at(record,hits.hits[1].u,hits.hits[1].v);
    require(a>=0&&b>=0&&a!=b,"Two roots were assigned to the same branch.");
    require(FoldedPatchGeometry::branch_at(record,.5,.5)==-1,"Exact fold point was silently assigned to a branch.");
}
void witness_test(const char* csv)
{
    const auto a=make_witness();const auto p=prepare(a);
    require(p.records.size()==12,"Expected twelve actual problem patches.");
    for(const auto& r:p.records)check_record(r);
    std::ofstream out;
    if(csv){out.open(csv);if(!out)throw std::runtime_error("Cannot create replay CSV.");
        out<<std::setprecision(17)<<"direction_id,old_flags,new_flags,folded_hits,focal_s,focal_p,incoherent_s,incoherent_p\n";}
    unsigned rescued=0,unchanged=0;
    for(unsigned i=0;i<a.directions.size();++i)
    {
        FocalOpticalResult old_f{},new_f{};
        const auto old=PatchOpticalEvaluator::evaluate_direction(a.view(),i,&a.focal,&old_f);
        const auto now=PatchOpticalEvaluator::evaluate_direction(a.view(),i,&a.focal,&new_f,p.view());
        require(now.known_hits_complete()&&new_f.valid(),"Captured direction was not fully evaluated.");
        require(std::isfinite(new_f.intensity_s+new_f.intensity_p),"Captured focal output remains nonfinite.");
        if(old.flags)
        {
            ++rescued;require(now.folded_evaluated_hits()==1,"Missing recovered contribution.");
            require(now.evaluated_hits==old.evaluated_hits+1,"Recovered field was not included in the sum.");
            require(now.regular_partial_incoherent_s+now.regular_partial_incoherent_p
                  >old.regular_partial_incoherent_s+old.regular_partial_incoherent_p,"No recovered incoherent energy.");
        }
        else
        {
            ++unchanged;
            // Structs have explicit initialized members and zeroed padding on
            // these targets; compare actual numeric/state members, not padding.
            require(now.regular_partial_incoherent_s==old.regular_partial_incoherent_s
                &&now.regular_partial_incoherent_p==old.regular_partial_incoherent_p
                &&now.regular_partial_path_s==old.regular_partial_path_s
                &&now.regular_partial_path_p==old.regular_partial_path_p
                &&new_f.intensity_s==old_f.intensity_s&&new_f.intensity_p==old_f.intensity_p,
                "Regular-only direction changed.");
        }
        if(out)out<<a.original_direction_ids[i]<<','<<old.flags<<','<<now.flags<<','<<now.folded_evaluated_hits()
            <<','<<new_f.intensity_s<<','<<new_f.intensity_p
            <<','<<now.regular_partial_incoherent_s<<','<<now.regular_partial_incoherent_p<<'\n';
    }
    require(rescued==12&&unchanged==24,"Unexpected recovery counts.");
    // Forced quadrature failure must not be marked successful or filled with 0.
    auto one=a.patches[p.records[0].compact_index];FoldedPatchConfig tiny{};tiny.maximum_panels=1;
    const auto failed=FoldedPatchBuilder::prepare(a.vertices.data(),static_cast<unsigned>(a.vertices.size()),one,
        p.records[0].compact_index,a.incident,&a.focal,tiny);
    require(failed.flags==FoldedQuadratureUnconverged,"Quadrature limit did not propagate failure.");
    // A malformed summary with an unexplained refinement flag is not certified.
    auto malformed=a;malformed.summaries[0].flags|=PatchQueryRefinementHit;
    FocalOpticalResult f{};const auto invalid=PatchOpticalEvaluator::evaluate_direction(malformed.view(),0,&a.focal,&f,p.view());
    require(!invalid.known_hits_complete(),"Unexplained upstream pending flag was cleared.");
    std::cout<<"Actual snapshot: 12/12 missing contributions restored; 24/24 regular-only directions unchanged.\n";
    if(out){out.flush();if(!out)throw std::runtime_error("Replay CSV write failed.");}
}
}
int main(int argc,char** argv)
{
    try
    {
        if(argc>2)throw std::invalid_argument("Usage: folded_cpu_test [replay.csv]");
        component_tests();double_root_test();witness_test(argc==2?argv[1]:nullptr);
        std::cout<<folded_tests::checks<<" folded-patch checks passed (CPU).\n";return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
