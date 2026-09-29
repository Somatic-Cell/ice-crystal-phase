#include "patch_validation.hpp"

#include <rainbow/raindrop_settings.hpp>
#include <cstdlib>
#include <iostream>
#include <limits>

int main()
{
    using namespace rainbow;
    using namespace rainbow::tests;
    try
    {
        auto geometry=square_patch();
        require_patch(geometry.certified_orientation()==1,"Square orientation is not positive.");
        const double exact=2.0*std::acos(-1.0)/3.0; // One cube face occupies one sixth of the sphere.
        require_patch(std::abs(geometry.signed_solid_angle()-exact)<3e-6,"Analytic solid angle failed.");
        require_patch(std::abs(integrated_patch_angle(geometry)-exact)<5e-5,"Independent integral failed.");
        std::swap(geometry.corners[0],geometry.corners[1]);
        std::swap(geometry.corners[2],geometry.corners[3]);
        require_patch(geometry.certified_orientation()==-1,"Negative orientation was discarded.");
        require_patch(std::abs(geometry.signed_solid_angle()+exact)<3e-6,"Signed area failed.");
        std::swap(geometry.corners[1],geometry.corners[3]);
        require_patch(geometry.certified_orientation()==0,"Fold was not flagged.");
        geometry=square_patch();
        for(auto& corner:geometry.corners) corner={1,0,0};
        require_patch(geometry.certified_orientation()==0,"Degenerate patch was certified.");

        RaindropSettings settings;settings.grid_width=2;settings.grid_height=2;
        PatchBuildLayout layout{};
        require_patch(PatchBuildLayout::try_make(settings.make_config(),layout),"Layout construction failed.");
        auto vertices=synthetic_patch_vertices();
        const PatchCellStatus expected[]={PatchCellStatus::RegularPositive,PatchCellStatus::NeedsRefinement,
                                         PatchCellStatus::MissingCorners,PatchCellStatus::NoOutgoingCorners};
        std::vector<PatchCellStatus> statuses(4);
        std::vector<OutgoingPatch> records;
        std::vector<PatchAabb> bounds;
        for(unsigned id=0;id<4;++id)
        {
            const auto result=PatchConstruction::make(vertices.data(),layout,id);
            require_patch(result.patch.status==expected[id],"Synthetic classification failed.");
            statuses[id]=result.patch.status;
            if(has_patch_geometry(statuses[id])) {records.push_back(result.patch);bounds.push_back(result.aabb);}
        }
        verify_patch_buffer(vertices,layout,statuses,records,bounds);
        // Distinct failure causes must survive neighboring absent corners.
        vertices[9].status=VertexStatus::PhaseOverflow;
        require_patch(PatchConstruction::make(vertices.data(),layout,2).patch.status==PatchCellStatus::InvalidVertex,
                      "Numerical failure was hidden by missing corner.");
        vertices[0].direction_drop.x=std::numeric_limits<float>::quiet_NaN();
        require_patch(PatchConstruction::make(vertices.data(),layout,0).patch.status==PatchCellStatus::InvalidGeometry,
                      "NaN direction was accepted.");

        // Independent analytic row-major mapping on a nonsquare grid.
        settings.grid_width=3;settings.grid_height=4;
        require_patch(PatchBuildLayout::try_make(settings.make_config(),layout),"Nonsquare layout failed.");
        require_patch(layout.vertices_per_path==12 && layout.cells_per_path==6 && layout.cell_count==24,
                      "Nonsquare counts failed.");
        std::uint32_t ids[4];layout.corner_indices(23,ids);
        require_patch(ids[0]==43 && ids[1]==44 && ids[2]==46 && ids[3]==47,"Corner IDs failed.");
        auto bad=settings.make_config();bad.grid_width=std::numeric_limits<std::uint32_t>::max();
        require_patch(!PatchBuildLayout::try_make(bad,layout),"Grid overflow was not rejected.");

        // Longitudes need no seam handling: geometry remains Cartesian.
        geometry=square_patch();
        for(auto& p:geometry.corners) p={-p.z,p.y,p.x};
        require_patch(geometry.certified_orientation()==1,"Longitude seam changed orientation.");
        require_patch(std::abs(geometry.signed_solid_angle()-exact)<3e-6,"Longitude seam changed area.");
        std::cout<<"Patch CPU geometry: passed (orientation, area, bounds, IDs, invalid/boundary cases).\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return EXIT_FAILURE;}
}
