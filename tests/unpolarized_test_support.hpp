#pragma once

#include <rainbow/raindrop_paths.hpp>
#include <rainbow/raindrop_settings.hpp>
#include <rainbow/patch_optical_evaluator.hpp>
#include <rainbow/folded_patch_builder.hpp>
#include <rainbow/query_direction_grid.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace unpolarized_tests
{
using namespace rainbow;
inline std::uint64_t checks = 0;
inline void require(bool condition, const char* message)
{ ++checks; if(!condition) throw std::runtime_error(message); }
inline bool near(double a, double b, double relative=5e-6, double absolute=1e-12)
{ return std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=absolute+relative*std::max(std::abs(a),std::abs(b)); }
inline bool near(Field32 a, Field32 b, double r=3e-6)
{ return near(a.x.real,b.x.real,r)&&near(a.x.imag,b.x.imag,r)&&near(a.y.real,b.y.real,r)&&near(a.y.imag,b.y.imag,r); }
inline bool near(OpticalField64 a, OpticalField64 b, double r=1e-5)
{ return near(a.s_real,b.s_real,r)&&near(a.s_imag,b.s_imag,r)&&near(a.p_real,b.p_real,r)&&near(a.p_imag,b.p_imag,r); }
inline bool same_geometry(const OutgoingVertex& a, const OutgoingVertex& b)
{
    return a.status==b.status && a.diagnostics==b.diagnostics
        && std::memcmp(&a.position_drop,&b.position_drop,sizeof(Vec3))==0
        && std::memcmp(&a.direction_drop,&b.direction_drop,sizeof(Vec3))==0
        && std::memcmp(&a.basis_x,&b.basis_x,sizeof(Vec3))==0
        && std::memcmp(&a.optical_cycles,&b.optical_cycles,sizeof(PhaseCycles))==0;
}
struct CountedIntersector
{
    RaindropIntersector solver;
    std::uint64_t calls=0;
    BoundaryHit intersect(const BoundaryRay& ray) {++calls;return solver.intersect(ray);}
};
struct Traced
{
    RaindropTraceConfig config;
    std::vector<OutgoingVertex> vertices;
    std::vector<Field32> second;
    std::uint64_t calls=0;
};
inline RaindropSettings settings(float radius, unsigned grid, bool sphere=false)
{
    RaindropSettings s;
    s.radius_mm=radius;s.grid_width=s.grid_height=grid;s.force_sphere=sphere;
    constexpr double a=20.0*QueryDirectionGrid::pi/180.0;
    s.incident_direction={float(std::cos(a)),float(-std::sin(a)),0};
    return s;
}
inline Traced trace(const RaindropSettings& s, bool paired)
{
    Traced a;a.config=s.make_config();const auto n=a.config.vertex_count();
    a.vertices.resize(std::size_t(n)*4);
    if(paired)a.second.resize(a.vertices.size());
    CountedIntersector intersection{{a.config.shape}};
    for(unsigned i=0;i<n;++i)
    {
        if(paired)RaindropPathTracer::trace_vertex_unpolarized(a.config,i,intersection,a.vertices.data(),a.second.data());
        else RaindropPathTracer::trace_vertex(a.config,i,intersection,a.vertices.data());
    }
    a.calls=intersection.calls;return a;
}
inline void check_trace(const Traced& pair, const Traced& x, const Traced& y)
{
    require(pair.calls==x.calls && pair.calls==y.calls,"Paired trace duplicated or changed geometric intersections.");
    require(pair.vertices.size()==x.vertices.size()&&pair.vertices.size()==y.vertices.size(),"Vertex counts differ.");
    for(std::size_t i=0;i<pair.vertices.size();++i)
    {
        require(same_geometry(pair.vertices[i],x.vertices[i])&&same_geometry(pair.vertices[i],y.vertices[i]),"Paired and separate traces changed geometry/status/path.");
        if(pair.vertices[i].status!=VertexStatus::Valid)continue;
        require(near(pair.vertices[i].field,x.vertices[i].field),"Input-column 0 differs from standalone x trace.");
        require(near(pair.second[i],y.vertices[i].field),"Input-column 1 differs from standalone y trace.");
    }
}
struct Geometry
{
    std::vector<OutgoingPatch> patches;
    std::vector<std::uint32_t> lookup;
    std::vector<FoldedPatchRecord> folded;
    FocalPhaseConfig focal{};
    FoldedPatchView view() const
    {return {lookup.data(),folded.data(),static_cast<unsigned>(lookup.size()),static_cast<unsigned>(folded.size())};}
};
inline Geometry build(const Traced& a)
{
    Geometry g;PatchBuildLayout layout{};
    require(PatchBuildLayout::try_make(a.config,layout),"Bad test layout.");
    g.focal.grid_width=a.config.grid_width;g.focal.grid_height=a.config.grid_height;
    g.focal.grid_half_extent=a.config.grid_half_extent;
    for(unsigned id=0;id<layout.cell_count;++id)
    {
        auto candidate=PatchConstruction::make(a.vertices.data(),layout,id);
        if(has_patch_geometry(candidate.patch.status))g.patches.push_back(candidate.patch);
    }
    g.lookup.assign(g.patches.size(),0xffffffffu);
    for(unsigned i=0;i<g.patches.size();++i)
    {
        if(g.patches[i].status!=PatchCellStatus::NeedsRefinement)continue;
        g.lookup[i]=static_cast<unsigned>(g.folded.size());
        g.folded.push_back(FoldedPatchBuilder::prepare(a.vertices.data(),static_cast<unsigned>(a.vertices.size()),
            g.patches[i],i,a.config.incident_direction,&g.focal));
    }
    return g;
}
struct Query
{
    Vec3 direction{};
    std::vector<PatchQueryHit> hits;
    PatchQuerySummary summary{};
    std::uint64_t offsets[2]{};
    PatchOpticsParams params(const Traced& a,const Geometry& g,bool dual) const
    {
        PatchOpticsParams p{};
        p.vertices=a.vertices.data();p.vertex_count=static_cast<unsigned>(a.vertices.size());
        p.patches=g.patches.data();p.patch_count=static_cast<unsigned>(g.patches.size());
        p.directions=&direction;p.direction_count=1;p.offsets=offsets;p.hits=hits.data();
        p.summaries=&summary;p.hit_storage_count=hits.size();
        p.incident_direction=a.config.incident_direction;p.incident_basis_x=a.config.incident_basis_x;
        if(dual){p.input_polarization=IncidentPolarization::Unpolarized;p.second_input_fields=a.second.data();p.second_input_count=p.vertex_count;}
        return p;
    }
};
inline Query query(const Traced& a,const Geometry& g,Vec3 direction)
{
    Query q;q.direction=direction;PatchHitCollector counter{};
    for(unsigned i=0;i<g.patches.size();++i)
    {
        BilinearPatchGeometry geometry{};
        for(unsigned k=0;k<4;++k)geometry.corners[k]=a.vertices[g.patches[i].vertex_indices[k]].direction_drop;
        // The same local primitive test is used twice just as the production count/fill passes.
        counter.consume(g.patches[i],i,BilinearPatchIntersector::intersect(geometry,direction));
    }
    q.hits.resize(counter.count);PatchHitCollector collector{};
    collector.storage=q.hits.data();collector.capacity=counter.count;collector.writing=true;
    for(unsigned i=0;i<g.patches.size();++i)
    {
        BilinearPatchGeometry geometry{};
        for(unsigned k=0;k<4;++k)geometry.corners[k]=a.vertices[g.patches[i].vertex_indices[k]].direction_drop;
        collector.consume(g.patches[i],i,BilinearPatchIntersector::intersect(geometry,direction));
    }
    q.summary=collector.finish(g.patches.data(),static_cast<unsigned>(g.patches.size()));
    q.offsets[1]=q.hits.size();return q;
}
inline void check_average(const PatchOpticalResult& p,const FocalOpticalResult& f,
                          const PatchOpticalResult& x,const FocalOpticalResult& fx,
                          const PatchOpticalResult& y,const FocalOpticalResult& fy,double tol=1e-5)
{
    require(p.flags==x.flags&&p.flags==y.flags,"Input averaging changed optical flags.");
    require(f.flags==fx.flags&&f.flags==fy.flags,"Input averaging changed focal flags.");
    require(p.hit_count==x.hit_count&&p.evaluated_hits==x.evaluated_hits&&p.rejected_hits==x.rejected_hits,"Hit counters were duplicated.");
    require(p.folded_evaluated_hits()==x.folded_evaluated_hits(),"Folded-hit count changed.");
    if((p.flags&patch_optical_error_mask)==0u)
    {
        require(near(p.regular_partial_path_field,x.regular_partial_path_field,tol),"Output first response column mismatch.");
        require(near(p.regular_partial_path_field_second,y.regular_partial_path_field,tol),"Output second response column mismatch.");
        require(near(p.regular_partial_incoherent_s,.5*(x.regular_partial_incoherent_s+y.regular_partial_incoherent_s),tol),"Unpolarized incoherent s mismatch.");
        require(near(p.regular_partial_incoherent_p,.5*(x.regular_partial_incoherent_p+y.regular_partial_incoherent_p),tol),"Unpolarized incoherent p mismatch.");
        require(near(p.regular_partial_path_s,.5*(x.regular_partial_path_s+y.regular_partial_path_s),tol),"Unpolarized coherent s mismatch.");
        require(near(p.regular_partial_path_p,.5*(x.regular_partial_path_p+y.regular_partial_path_p),tol),"Unpolarized coherent p mismatch.");
    }
    if(f.valid())
    {
        require(near(f.intensity_s,.5*(fx.intensity_s+fy.intensity_s),tol)&&near(f.intensity_p,.5*(fx.intensity_p+fy.intensity_p),tol),"Unpolarized focal mismatch.");
        require(near(f.field,fx.field,tol)&&near(f.field_second,fy.field,tol),"Focal response columns mismatch.");
        for(unsigned k=0;k<4;++k)
        {
            require(f.family_hits[k]==fx.family_hits[k],"Family counts doubled.");
            require(near(f.family_incoherent[k],.5*(fx.family_incoherent[k]+fy.family_incoherent[k]),tol),"Family intensity is not averaged.");
        }
    }
}
}
