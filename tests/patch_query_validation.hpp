#pragma once

#include <rainbow/patch_query_data.hpp>
#include <rainbow/query_direction_grid.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace rainbow::tests
{
inline void require_query(const bool ok,const char* message)
{if(!ok)throw std::runtime_error(message);}

// 独立した参照: long double の3次元 scalar triple products で係数を作る．
// 本番の2次元 projection / difference_of_products / Newton は呼ばない．
// MSVC の long double は double と同精度だが，FP32 本番演算とは独立である．
struct ReferenceBilinearQuery
{
    struct V
    {
        long double x,y,z;
        V operator+(V b)const{return{x+b.x,y+b.y,z+b.z};}
        V operator-(V b)const{return{x-b.x,y-b.y,z-b.z};}
        V operator*(long double a)const{return{x*a,y*a,z*a};}
        V cross(V b)const{return{y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x};}
        long double dot(V b)const{return x*b.x+y*b.y+z*b.z;}
    };
    struct Hit{long double u,v,t;};
    static V convert(Vec3 a){return{a.x,a.y,a.z};}
    static std::vector<Hit> solve(const BilinearPatchGeometry& g,Vec3 direction)
    {
        const auto w=convert(direction),p0=convert(g.corners[0]),p1=convert(g.corners[1]);
        const auto p2=convert(g.corners[2]),p3=convert(g.corners[3]);
        const auto e=p1-p0,f=p2-p0,h=p3-p2-e;
        const long double a=w.dot(e.cross(h));
        const long double b=w.dot(p0.cross(h)+e.cross(f));
        const long double c=w.dot(p0.cross(f));
        std::vector<long double> roots;
        if(a==0){if(b!=0)roots.push_back(-c/b);}
        else
        {
            const long double discriminant=b*b-4*a*c;
            if(discriminant==0)roots.push_back(-b/(2*a));
            if(discriminant>0)
            {
                const auto q=-0.5L*(b+std::copysign(std::sqrt(discriminant),b));
                roots={q/a,c/q};
            }
        }
        std::vector<Hit> result;
        for(auto u:roots)
        {
            if(u<0 || u>1)continue;
            // Intersect the ray with the v-isoline at this u (different recovery from production).
            const auto o=p0+e*u,d=f+h*u,n=w.cross(d);
            const auto den=n.dot(n);
            if(den==0)continue;
            const auto v=o.cross(w).dot(n)/den;
            const auto t=o.cross(d).dot(n)/den;
            if(v>=0 && v<=1 && t>0)result.push_back({u,v,t});
        }
        std::sort(result.begin(),result.end(),[](auto a0,auto b0){return a0.u<b0.u;});
        return result;
    }
};

inline BilinearPatchGeometry front_patch()
{
    BilinearPatchGeometry p{{{-0.2f,-0.2f,1.0f},{0.2f,-0.2f,1.0f},
                            {-0.2f,0.2f,1.0f},{0.2f,0.2f,1.0f}}};
    for(auto& v:p.corners)v=v.normalized();
    return p;
}
inline BilinearPatchGeometry folded_patch()
{
    BilinearPatchGeometry p{};
    const float xs[4]={-0.5f,0,0,0.5f},ys[4]={-0.08f,-0.08f,-0.08f,0.42f};
    for(unsigned i=0;i<4;++i)p.corners[i]={xs[i],ys[i],std::sqrt(1-xs[i]*xs[i]-ys[i]*ys[i])};
    return p;
}
inline BilinearPatchGeometry two_depth_patch()
{
    // 1本の方向が同じ非平面パッチを2回通る．近い方だけ返す実装を検出する．
    return {{{0.384129852f,-0.294093579f,0.875187576f},
             {-0.642671168f,0.357146144f,0.677805662f},
             {-0.255662054f,0.0962326676f,0.961964726f},
             {0.919800222f,-0.187220544f,0.344841957f}}};
}
inline std::vector<OutgoingVertex> query_test_vertices()
{
    std::array<BilinearPatchGeometry,4> patches{front_patch(),front_patch(),two_depth_patch(),front_patch()};
    std::swap(patches[1].corners[0],patches[1].corners[1]);
    std::swap(patches[1].corners[2],patches[1].corners[3]);
    for(auto& p:patches[3].corners)p.z=-p.z;
    std::vector<OutgoingVertex> vertices(16);
    for(unsigned family=0;family<4;++family)for(unsigned i=0;i<4;++i)
    {
        auto& v=vertices[family*4+i];
        v.direction_drop=patches[family].corners[i];
        v.basis_x=TransverseFrame::from_direction(v.direction_drop).e0;
        v.field.x={1,0};v.status=VertexStatus::Valid;
    }
    return vertices;
}

inline void verify_hit(const BilinearPatchGeometry& p,Vec3 direction,const BilinearRayHit& h)
{
    require_query(h.u>=0 && h.u<=1 && h.v>=0 && h.v<=1 && h.t>0,"Invalid (u,v,t).");
    const auto e=p.evaluate(h.u,h.v)-direction*h.t;
    require_query(e.length()<2e-6f,"Ray/bilinear reconstruction residual too large.");
    const auto box=p.bounds();const auto r=direction*h.t;
    require_query(r.x>=box.min_x&&r.x<=box.max_x&&r.y>=box.min_y&&r.y<=box.max_y
        &&r.z>=box.min_z&&r.z<=box.max_z,"Reported hit lies outside the built AABB.");
}

inline void test_bilinear_queries()
{
    const auto p=front_patch();
    auto r=BilinearPatchIntersector::intersect(p,{0,0,1});
    require_query(r.count==1 && r.flags==0,"Planar center intersection.");
    require_query(std::abs(r.hits[0].u-.5f)<1e-6f && std::abs(r.hits[0].v-.5f)<1e-6f,"Planar center uv.");
    verify_hit(p,{0,0,1},r.hits[0]);
    require_query(BilinearPatchIntersector::intersect(p,{0,0,-1}).count==0,"Behind-ray rejection.");
    require_query(BilinearPatchIntersector::intersect(p,{1,0,0}).count==0,"Parallel outside rejection.");
    auto edge=BilinearPatchIntersector::intersect(p,p.corners[0]);
    require_query(edge.count==1 && (edge.hits[0].flags&BilinearBoundary),"Closed corner must be retained and flagged.");
    r=BilinearPatchIntersector::intersect(folded_patch(),{0,0,1});
    require_query(r.count==2 && !(r.flags&BilinearUnresolved),"Both folded-map roots are required.");
    require_query(std::abs(r.hits[0].u-.2f)<2e-6f && std::abs(r.hits[1].u-.8f)<2e-6f,"Folded roots.");
    for(unsigned i=0;i<2;++i)verify_hit(folded_patch(),{0,0,1},r.hits[i]);
    const auto two_depth=BilinearPatchIntersector::intersect(two_depth_patch(),{0,0,1});
    require_query(two_depth.count==2 && std::abs(two_depth.hits[0].t-two_depth.hits[1].t)>.05f,
        "Two distinct depths must both be collected.");
    BilinearPatchGeometry point{};for(auto& c:point.corners)c={0,0,1};
    require_query((BilinearPatchIntersector::intersect(point,{0,0,1}).flags&BilinearUnresolved)!=0,
        "Non-isolated solution must not be silently called a miss.");
    require_query((BilinearPatchIntersector::intersect(p,{0,0,0}).flags&BilinearInvalidInput)!=0,"Zero ray direction.");

    
    // 丸めなく表現できる人工的なパッチで，重解の処理を検査する．
    // 雨粒の形状や本番パッチを置き換えるものではない．
    //
    // P(u,v) = (u + v - 1, (u - 0.5)(v - 0.5), 1)
    // +z 方向との交差は u = v = 0.5, t = 1 の一つ．
    const BilinearPatchGeometry repeated_patch{{
        {-1.0f,  0.25f, 1.0f},
        { 0.0f, -0.25f, 1.0f},
        { 0.0f, -0.25f, 1.0f},
        { 1.0f,  0.25f, 1.0f}
    }};
    
    const auto repeated_result =
        BilinearPatchIntersector::intersect(
            repeated_patch, Vec3{0.0f, 0.0f, 1.0f});
        
    require_query(
        repeated_result.count == 1
            && (repeated_result.flags & BilinearUnresolved) == 0u,
        "Repeated root must be returned once.");
    
    const auto& repeated_hit = repeated_result.hits[0];
    
    require_query(
        repeated_hit.u == 0.5f
            && repeated_hit.v == 0.5f
            && repeated_hit.t == 1.0f,
        "Repeated root coordinates.");
    
    require_query(
        (repeated_hit.flags & BilinearSingular) != 0u
            && (repeated_hit.flags & BilinearUsedFp64) != 0u,
        "Repeated root must retain singular and FP64 flags.");
    // Independent reference for non-edge queries, both arbitrary small patches
    // and broad folded patches. Input ray rounding is part of the exact test input.
    std::mt19937 engine(0x119f70u);
    std::uniform_real_distribution<float> uniform(-1,1),inside(.002f,.998f);
    for(unsigned test=0;test<12000;++test)
    {
        const float scale=std::pow(10.0f,-1.0f-float(test%4));
        const float x=uniform(engine)*.4f,y=uniform(engine)*.4f;
        BilinearPatchGeometry g{{{x-scale,y-scale,1},{x+scale,y-scale,1},
            {x-scale,y+scale,1},{x+scale,y+scale,1+scale*.2f*uniform(engine)}}};
        for(auto& c:g.corners)c=c.normalized();
        Vec3 w{};
        if(test%3)w=g.evaluate(inside(engine),inside(engine)).normalized();
        else w=Vec3{uniform(engine),uniform(engine),1}.normalized();
        const auto got=BilinearPatchIntersector::intersect(g,w);
        const auto reference=ReferenceBilinearQuery::solve(g,w);
        if(got.count!=reference.size() || (got.flags&BilinearUnresolved))
        {std::ostringstream m;m<<"Reference root count mismatch, case "<<test;throw std::runtime_error(m.str());}
        for(unsigned i=0;i<got.count;++i)
        {
            verify_hit(g,w,got.hits[i]);
            require_query(std::abs((long double)got.hits[i].u-reference[i].u)<3e-5L
                &&std::abs((long double)got.hits[i].v-reference[i].v)<3e-5L,"Independent reference uv mismatch.");
        }
    }
}
inline void test_query_collector()
{
    OutgoingPatch patches[3]{};
    patches[0].patch_id=10;patches[1].patch_id=8;patches[2].patch_id=9;
    patches[0].status=PatchCellStatus::RegularPositive;
    patches[1].status=PatchCellStatus::RegularNegative;
    patches[2].status=PatchCellStatus::NeedsRefinement;
    const auto single=BilinearPatchIntersector::intersect(front_patch(),{0,0,1});
    const auto two=BilinearPatchIntersector::intersect(folded_patch(),{0,0,1});
    PatchHitCollector counter{};
    counter.consume(patches[0],0,single);counter.consume(patches[2],2,two);
    counter.consume(patches[1],1,single);
    require_query(counter.count==4,"Count pass must include both roots.");
    PatchQueryHit storage[5]{};
    PatchHitCollector writer{storage,5,0,0,true};
    writer.consume(patches[2],2,two);writer.consume(patches[0],0,single);
    writer.consume(patches[1],1,single);writer.consume(patches[0],0,single); // duplicate primitive
    const auto result=writer.finish(patches,3);
    require_query(result.hit_count==4 && result.regular_hits==2 && result.refinement_hits==2,"Collector classification.");
    require_query((result.flags&PatchQueryDuplicateReport)!=0,"Duplicate report flag.");
    require_query(storage[0].patch_id==8 && storage[1].patch_id==9 && storage[2].root_index==1
        &&storage[3].patch_id==10,"Deterministic ordering.");
    PatchQueryHit guard[2]{};guard[1].patch_id=0xdeadbeefu;
    PatchHitCollector short_buffer{guard,1,0,0,true};
    short_buffer.consume(patches[2],2,two);
    require_query((short_buffer.finish(patches,3).flags&PatchQueryCountMismatch)!=0
        &&guard[1].patch_id==0xdeadbeefu,"Capacity guard must not overwrite memory.");
}
inline void test_direction_grid()
{
    RaindropTraceConfig c{};
    c.incident_direction=Vec3{1,-.25f,0}.normalized();
    c.incident_basis_x=TransverseFrame::from_direction(c.incident_direction).e0;
    const auto g=QueryDirectionGrid::make(c,71,113);const auto d=g.directions();
    double total=0;
    for(std::uint32_t row=0;row<g.theta_count;++row)total+=g.solid_angle(row)*g.phi_count;
    require_query(std::abs(total-4*QueryDirectionGrid::pi)<1e-12,"Longitude-latitude solid angle weights.");
    for(std::size_t i=0;i<d.size();++i)
        require_query(std::abs(d[i].dot(c.incident_direction)-std::cos(g.theta(std::uint32_t(i/g.phi_count))))<3e-7,
            "Query angular convention.");
}
} // namespace rainbow::tests
