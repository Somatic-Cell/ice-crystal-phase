#pragma once
#include <rainbow/folded_patch_builder.hpp>
#include <rainbow/patch_optical_evaluator.hpp>
#include "data/folded_witness.hpp"
#include <cmath>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace folded_tests
{
inline unsigned checks=0;
inline void require(bool c,const char* message)
{++checks;if(!c)throw std::runtime_error(message);}
inline bool near(double a,double b,double r=3e-7,double absolute=1e-12)
{return std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=absolute+r*std::max(std::abs(a),std::abs(b));}
struct Prepared
{
    std::vector<std::uint32_t> lookup;
    std::vector<FoldedPatchRecord> records;
    FoldedPatchView view() const {return {lookup.data(),records.data(),static_cast<unsigned>(lookup.size()),static_cast<unsigned>(records.size())};}
};
inline Prepared prepare(const Witness& a)
{
    Prepared p;p.lookup.assign(a.patches.size(),0xffffffffu);
    for(auto index:a.used_patches)
    {
        const auto& patch=a.patches[index];
        if(patch.status!=PatchCellStatus::NeedsRefinement)continue;
        p.lookup[index]=static_cast<unsigned>(p.records.size());
        p.records.push_back(FoldedPatchBuilder::prepare(a.vertices.data(),static_cast<unsigned>(a.vertices.size()),
            patch,index,a.incident,&a.focal));
    }
    return p;
}
inline void check_record(const FoldedPatchRecord& r)
{
    require(r.flags==FoldedReady,"Actual witness branch preparation failed.");
    require(r.branch_count==2,"Actual witness did not produce two branches.");
    double area=0,energy=0;
    for(unsigned j=0;j<r.branch_count;++j)
    {
        const auto& b=r.branches[j];area+=b.area_fraction;
        require(b.area_fraction>0&&b.solid_angle_sr>0,"Nonpositive branch measure.");
        require(b.focal_flags==0,"Actual witness focal representative unresolved.");
        require(FoldedPatchGeometry::branch_at(r,b.representative_u,b.representative_v)==static_cast<int>(j),"Representative not inside its branch.");
        // Finite uniform branch density integrates to its assigned input area.
        energy+=(b.area_fraction/b.solid_angle_sr)*b.solid_angle_sr;
    }
    require(near(area,1,1e-9,1e-12),"Input area duplicated or lost.");
    require(near(energy,1,1e-9,1e-12),"Constant-envelope energy allocation changed.");
    if(r.patch_id==54155)
    {
        require(near(r.branches[0].area_fraction,0.868511199133002,1e-9),"Independent parameter area mismatch.");
        require(near(r.branches[0].solid_angle_sr,5.007818001209945e-5,1e-9,1e-15),"Independent negative-branch solid angle mismatch.");
        require(near(r.branches[1].solid_angle_sr,2.175586525822955e-6,1e-9,1e-15),"Independent positive-branch solid angle mismatch.");
    }
}
inline void check_result_pair(const PatchOpticalResult& a,const FocalOpticalResult& f,
                             const PatchOpticalResult& reference,const FocalOpticalResult& ref_f)
{
    require(a.flags==reference.flags&&a.evaluated_hits==reference.evaluated_hits
       &&a.rejected_hits==reference.rejected_hits&&a.folded_evaluated_hits()==reference.folded_evaluated_hits(),"Optical state/count mismatch.");
    require(f.flags==ref_f.flags&&f.corrected_hits==ref_f.corrected_hits,"Focal state/count mismatch.");
    require(near(a.regular_partial_incoherent_s,reference.regular_partial_incoherent_s)
        &&near(a.regular_partial_incoherent_p,reference.regular_partial_incoherent_p),"Incoherent CPU/GPU mismatch.");
    require(near(a.regular_partial_path_s,reference.regular_partial_path_s,1e-5)
        &&near(a.regular_partial_path_p,reference.regular_partial_path_p,1e-5),"Path CPU/GPU mismatch.");
    require(near(f.intensity_s,ref_f.intensity_s,1e-5)&&near(f.intensity_p,ref_f.intensity_p,1e-5),"Focal CPU/GPU mismatch.");
}
}
