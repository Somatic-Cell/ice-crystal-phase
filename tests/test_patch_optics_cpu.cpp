#include "patch_optics_test_cases.hpp"
#include <cstdlib>
#include <iostream>

int main()
{
    try
    {
        using namespace rainbow;
        const auto input=tests::make_optics_fixture();
        const auto params=input.view();
        std::vector<PatchOpticalResult> output(input.directions.size());
        for(std::uint32_t i=0;i<params.direction_count;++i)
            output[i]=PatchOpticalEvaluator::evaluate_direction(params,i);
        tests::verify_optics_fixture(input,output);
        // Guard every indirect range before dereferencing hit/vertex pointers.
        auto bad=input;
        bad.offsets[1]=bad.hits.size()+1u;
        auto r=PatchOpticalEvaluator::evaluate_direction(bad.view(),0);
        tests::require_optics((r.flags & PatchOpticalInvalidInput)!=0u,"Out-of-range segment accepted.");
        bad=input;bad.hits[0].compact_index=0xffffffffu;
        r=PatchOpticalEvaluator::evaluate_direction(bad.view(),0);
        tests::require_optics((r.flags & PatchOpticalInvalidInput)!=0u,"Out-of-range patch accepted.");
        bad=input;bad.patches[0].vertex_indices[0]=0xffffffffu;
        r=PatchOpticalEvaluator::evaluate_direction(bad.view(),0);
        tests::require_optics((r.flags & PatchOpticalInvalidInput)!=0u,"Out-of-range vertex accepted.");
        bad=input;
        const auto start=bad.offsets[1];
        bad.hits[start+1u]=bad.hits[start];
        r=PatchOpticalEvaluator::evaluate_direction(bad.view(),1);
        tests::require_optics((r.flags & PatchOpticalInvalidInput)!=0u,"Duplicate known hit summed twice.");
        // Shared source-frame gauge changes must not create optical interference.
        auto gauge=input;
        constexpr float alpha=.731f;
        const float c=std::cos(alpha),s=std::sin(alpha);
        for(auto& vertex:gauge.vertices)
        {
            if(vertex.basis_x.z!=0 || vertex.direction_drop.z!=1.0f
               || !std::isfinite(vertex.field.x.real)) continue;
            const auto e0=vertex.basis_x;
            const auto e1=vertex.direction_drop.cross(e0).normalized();
            vertex.basis_x=e0*c+e1*s;
            vertex.field=vertex.field.rotated_basis(c,s);
        }
        // Compare only evaluable inputs, excluding the huge-cancellation stress
        // test where perturbing stored float components changes the input waves.
        for(std::uint32_t i=0;i<params.direction_count;++i)
            if(i!=8u && output[i].flags==0u)
                tests::compare_optics_results(PatchOpticalEvaluator::evaluate_direction(gauge.view(),i),output[i]);
        std::cout<<"Patch optics CPU: "<<output.size()<<" cases, dyadic reference, bounds, partial/NaN, gauge passed.\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return EXIT_FAILURE;}
}
