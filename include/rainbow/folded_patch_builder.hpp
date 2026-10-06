#pragma once

#include <rainbow/folded_patch_geometry.hpp>
#include <rainbow/focal_phase.hpp>

namespace rainbow
{
// Finite branch-density extension. Adaptive quadrature subdivides ONLY the
// integration rectangle, never the original patch or the incident ray grid.
struct FoldedPatchBuilder
{
    struct Integral
    {
        double area=0, moment_u=0, moment_v=0, omega=0;
        HOST_DEVICE void add(const Integral& b) noexcept
        {area+=b.area;moment_u+=b.moment_u;moment_v+=b.moment_v;omega+=b.omega;}
    };
    struct Panel {double u0,u1,z0,z1;unsigned depth;};

    [[nodiscard]] HOST_DEVICE static bool config_valid(const FoldedPatchConfig& c) noexcept
    {
        return c.relative_tolerance>=1e-12 && c.relative_tolerance<=1e-3
            && c.maximum_panels>=1 && c.maximum_panels<=65536
            && c.maximum_depth>=1 && c.maximum_depth<=24;
    }
    [[nodiscard]] HOST_DEVICE static Integral quadrature(const FoldedPatchGeometry& g,
        int sign,const Panel& p,bool high) noexcept
    {
        const double x8[8] = {
            -0.96028985649753618,
            -0.79666647741362673,
            -0.52553240991632899,
            -0.18343464249564978,
            0.18343464249564978,
            0.52553240991632899,
            0.79666647741362673,
            0.96028985649753618};
        const double w8[8] = {
            0.10122853629037706,
            0.22238103445337443,
            0.31370664587788688,
            0.36268378337836166,
            0.36268378337836166,
            0.31370664587788688,
            0.22238103445337443,
            0.10122853629037706};
        const double x16[16] = {
            -0.98940093499164994,
            -0.9445750230732326,
            -0.86563120238783176,
            -0.755404408355003,
            -0.61787624440264377,
            -0.45801677765722737,
            -0.28160355077925892,
            -0.095012509837637441,
            0.095012509837637441,
            0.28160355077925892,
            0.45801677765722737,
            0.61787624440264377,
            0.755404408355003,
            0.86563120238783176,
            0.9445750230732326,
            0.98940093499164994};
        const double w16[16] = {
            0.027152459411754176,
            0.062253523938647456,
            0.095158511682492605,
            0.12462897125553407,
            0.14959598881657671,
            0.16915651939500265,
            0.18260341504492364,
            0.18945061045506864,
            0.18945061045506864,
            0.18260341504492364,
            0.16915651939500265,
            0.14959598881657671,
            0.12462897125553407,
            0.095158511682492605,
            0.062253523938647456,
            0.027152459411754176};

        const unsigned n=high?16u:8u;
        const double* xs=high?x16:x8;const double* ws=high?w16:w8;
        const double uh=(p.u1-p.u0)*0.5, uc=(p.u0+p.u1)*0.5;
        const double zh=(p.z1-p.z0)*0.5, zc=(p.z0+p.z1)*0.5;
        Integral total{};
        for(unsigned i=0;i<n;++i)
        {
            const double u=::fma(uh,xs[i],uc);double vl=0,vh=0;
            if(!FoldedPatchGeometry::section(g.j,sign,u,vl,vh))continue;
            const double width=vh-vl;
            for(unsigned k=0;k<n;++k)
            {
                const double z=::fma(zh,xs[k],zc),v=::fma(z,width,vl);
                const double w=uh*zh*ws[i]*ws[k]*width;
                const auto point=g.point(u,v);
                const double q=point.dot(point);
                const double J=::fabs(FoldedPatchGeometry::value(g.j,u,v))*g.scale;
                total.area+=w;total.moment_u+=w*u;total.moment_v+=w*v;
                total.omega+=w*J/(q*::sqrt(q));
            }
        }
        return total;
    }
    [[nodiscard]] HOST_DEVICE static bool integrate(const FoldedPatchGeometry& g,FoldedBranch& b,
        const FoldedPatchConfig& c,FoldedPatchRecord& record) noexcept
    {
        double cuts[4]{};const unsigned n=g.cuts(cuts);
        Integral sum{};
        // Splitting at changes of section topology keeps the root curve smooth
        // inside every starting panel, including arbitrarily small components.
        for(unsigned k=0;k+1<n;++k)
        {
            const double lo=::fmax(cuts[k],b.u0),hi=::fmin(cuts[k+1],b.u1);
            if(!(hi>lo))continue;
            Panel stack[26]{};unsigned size=1;
            stack[0]={lo,hi,0,1,0};
            while(size)
            {
                const Panel p=stack[--size];
                if(++record.panels>c.maximum_panels)return false;
                const Integral low=quadrature(g,b.sign,p,false), high=quadrature(g,b.sign,p,true);
                if(!(high.area>0 && high.omega>0
                    && FoldedPatchGeometry::finite(high.area) && FoldedPatchGeometry::finite(high.omega)))return false;
                const double ea=::fabs(high.area-low.area)/high.area;
                const double es=::fabs(high.omega-low.omega)/high.omega;
                const double emu=::fabs(high.moment_u-low.moment_u)/high.area;
                const double emv=::fabs(high.moment_v-low.moment_v)/high.area;
                const double error=::fmax(::fmax(ea,es),::fmax(emu,emv));
                if(error<=c.relative_tolerance)
                {
                    record.largest_estimated_relative_error=::fmax(record.largest_estimated_relative_error,error);
                    sum.add(high);continue;
                }
                if(p.depth>=c.maximum_depth || size+2>26)return false;
                if((p.depth&1u)==0)
                {
                    const double m=(p.u0+p.u1)*0.5;
                    if(!(m>p.u0&&m<p.u1))return false;
                    stack[size++]={m,p.u1,p.z0,p.z1,p.depth+1};
                    stack[size++]={p.u0,m,p.z0,p.z1,p.depth+1};
                }
                else
                {
                    const double m=(p.z0+p.z1)*0.5;
                    if(!(m>p.z0&&m<p.z1))return false;
                    stack[size++]={p.u0,p.u1,m,p.z1,p.depth+1};
                    stack[size++]={p.u0,p.u1,p.z0,m,p.depth+1};
                }
            }
        }
        if(!(sum.area>0 && sum.omega>0 && FoldedPatchGeometry::finite(sum.omega)))return false;
        b.area_fraction=sum.area;b.solid_angle_sr=sum.omega;
        b.representative_u=sum.moment_u/sum.area;
        b.representative_v=sum.moment_v/sum.area;
        // A component's centroid need not be inside a nonconvex sign region.
        // If needed, use the midpoint of the vertical section at the u-centroid.
        double vl=0,vh=0;
        if(!(b.representative_u>b.u0&&b.representative_u<b.u1)
            || !FoldedPatchGeometry::section(g.j,b.sign,b.representative_u,vl,vh))return false;
        if(!(b.representative_v>vl&&b.representative_v<vh))b.representative_v=(vl+vh)*0.5;
        return b.representative_v>vl && b.representative_v<vh;
    }
    [[nodiscard]] HOST_DEVICE static FoldedPatchRecord geometry(
        const Vec3* corners,std::uint32_t patch_id,std::uint32_t compact_index,
        const FoldedPatchConfig& config={}) noexcept
    {
        FoldedPatchRecord r{};r.patch_id=patch_id;r.compact_index=compact_index;
        if(!corners || !config_valid(config))return r;
        FoldedPatchGeometry g{};r.flags=g.initialize(corners);
        if(r.flags!=FoldedReady)return r;
        for(unsigned i=0;i<4;++i)r.jacobian[i]=g.j[i];
        if(!g.components(r)){r.flags=FoldedTopologyUnresolved;return r;}
        double area=0;
        for(unsigned i=0;i<r.branch_count;++i)
        {
            if(!integrate(g,r.branches[i],config,r))
            {r.flags=FoldedQuadratureUnconverged;return r;}
            area+=r.branches[i].area_fraction;
        }
        // This is a consistency check, not a normalization. Do not rescale
        // branches to force conservation after inaccurate integration.
        if(!FoldedPatchGeometry::finite(area)||::fabs(area-1)>8*config.relative_tolerance)
            r.flags=FoldedInvalidMetric;
        return r;
    }
    [[nodiscard]] HOST_DEVICE static FoldedPatchRecord prepare(
        const OutgoingVertex* vertices,std::uint32_t vertex_count,
        const OutgoingPatch& patch,std::uint32_t compact_index,Vec3 incident,
        const FocalPhaseConfig* focal,const FoldedPatchConfig& config={}) noexcept
    {
        FoldedPatchRecord bad{};bad.patch_id=patch.patch_id;bad.compact_index=compact_index;
        if(!vertices || patch.status!=PatchCellStatus::NeedsRefinement
            || !(patch.incident_area_drop2>0) || !FoldedPatchGeometry::finite(patch.incident_area_drop2))return bad;
        Vec3 corners[4]{};
        for(unsigned i=0;i<4;++i)
        {
            const auto id=patch.vertex_indices[i];
            if(id>=vertex_count || vertices[id].status!=VertexStatus::Valid
                || !PatchConstruction::usable_vertex(vertices[id]))return bad;
            corners[i]=vertices[id].direction_drop;
        }
        auto r=geometry(corners,patch.patch_id,compact_index,config);
        if(r.flags!=FoldedReady)return r;
        if(focal)for(unsigned i=0;i<r.branch_count;++i)
        {
            auto& b=r.branches[i];FocalPhaseEstimate estimate{};
            b.focal_flags=FocalPhase::estimate_at(vertices,vertex_count,patch,incident,*focal,
                                                  b.representative_u,b.representative_v,estimate);
            if(b.focal_flags==FocalNone)
            {b.quarter_turns=estimate.quarter_turns;b.extra_quarter_turn=estimate.extra_quarter_turn;}
        }
        return r;
    }
};
} // namespace rainbow
