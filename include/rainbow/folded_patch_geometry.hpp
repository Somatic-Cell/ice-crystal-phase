#pragma once

#include <rainbow/folded_patch_data.hpp>
#include <rainbow/patch_geometry.hpp>
#include <cmath>

namespace rainbow
{
// Uses the ORIGINAL four chord vertices. Splitting integration domains does
// not tessellate, normalize, retrace, or create new optical samples.
struct FoldedPatchGeometry
{
    using D3 = Vec3T<double>;
    D3 corners[4]{};
    double j[4]{}, scale = 0;

    [[nodiscard]] HOST_DEVICE static bool finite(double x) noexcept
    { return x >= -0x1.fffffffffffffp1023 && x <= 0x1.fffffffffffffp1023; }
    [[nodiscard]] HOST_DEVICE static double mix(double a, double b, double t) noexcept
    { return ::fma(t, b-a, a); }
    [[nodiscard]] HOST_DEVICE static double value(const double* a, double u, double v) noexcept
    { return mix(mix(a[0],a[1],u),mix(a[2],a[3],u),v); }
    [[nodiscard]] HOST_DEVICE D3 point(double u,double v) const noexcept
    { return (corners[0]*(1-u)+corners[1]*u)*(1-v)+(corners[2]*(1-u)+corners[3]*u)*v; }
    [[nodiscard]] HOST_DEVICE static double determinant(D3 a,D3 b,D3 c) noexcept
    {
        const double yz=b.z*c.y, zx=b.x*c.z, xy=b.y*c.x;
        const double x=::fma(b.y,c.z,-yz)+::fma(-b.z,c.y,yz);
        const double y=::fma(b.z,c.x,-zx)+::fma(-b.x,c.z,zx);
        const double z=::fma(b.x,c.y,-xy)+::fma(-b.y,c.x,xy);
        return ::fma(a.x,x,::fma(a.y,y,a.z*z));
    }
    [[nodiscard]] HOST_DEVICE std::uint32_t initialize(const Vec3* input) noexcept
    {
        for(unsigned i=0;i<4;++i)
        {
            if(!input[i].is_finite()) return FoldedInvalidInput;
            corners[i]=input[i].cast<double>();
            const double q=corners[i].dot(corners[i]);
            if(!(q>0 && q<=4)) return FoldedInvalidInput;
        }
        const D3 axis=(corners[0]+corners[1])+(corners[2]+corners[3]);
        // Conservative FP64 roundoff guard for this stored-corner sum/dot.
        // Failure is NOT interpreted as proof of an origin crossing.
        for(unsigned i=0;i<4;++i)
        {
            double bound=0;
            for(unsigned k=0;k<4;++k)
                bound+=::fabs(corners[k].x*corners[i].x)+::fabs(corners[k].y*corners[i].y)
                      +::fabs(corners[k].z*corners[i].z);
            if(!(axis.dot(corners[i])>0x1p-46*bound)) return FoldedHemisphereUncertified;
        }
        const D3 eu=corners[1]-corners[0],ev=corners[2]-corners[0];
        j[0]=determinant(corners[0],eu,ev);
        j[1]=determinant(corners[1],eu,corners[3]-corners[1]);
        j[2]=determinant(corners[2],corners[3]-corners[2],ev);
        j[3]=determinant(corners[3],corners[3]-corners[2],corners[3]-corners[1]);
        scale=0;
        bool pos=false,neg=false;
        for(unsigned i=0;i<4;++i)
        {
            if(!finite(j[i])) return FoldedInvalidInput;
            scale=::fmax(scale,::fabs(j[i]));pos|=j[i]>0;neg|=j[i]<0;
        }
        if(!(pos&&neg&&scale>0)) return FoldedNoSignChange;
        for(unsigned i=0;i<4;++i) j[i]/=scale;
        return FoldedReady;
    }
    // Interval in v for sign*J(u,v)>0. Endpoints J=0 have measure zero.
    [[nodiscard]] HOST_DEVICE static bool section(const double* a,int sign,double u,
                                                  double& lo,double& hi) noexcept
    {
        const double b=double(sign)*mix(a[0],a[1],u);
        const double t=double(sign)*mix(a[2],a[3],u);
        if(b<=0 && t<=0) return false;
        lo=0;hi=1;
        if(b>0 && t>0) return true;
        // Opposite/nonpositive endpoint: this ratio is in [0,1] without a
        // subtractive cancellation. No epsilon removal of narrow components.
        const double root=::fabs(b)/(::fabs(b)+::fabs(t));
        if(b>0)hi=root;else lo=root;
        return hi>lo;
    }
    [[nodiscard]] HOST_DEVICE unsigned cuts(double* x) const noexcept
    {
        unsigned n=2;x[0]=0;x[1]=1;
        for(unsigned edge=0;edge<2;++edge)
        {
            const double a=j[2*edge],b=j[2*edge+1];
            if((a<0&&b>0)||(a>0&&b<0))
                x[n++]=::fabs(a)/(::fabs(a)+::fabs(b));
        }
        for(unsigned i=1;i<n;++i)
        {double t=x[i];unsigned k=i;while(k&&x[k-1]>t){x[k]=x[k-1];--k;}x[k]=t;}
        unsigned m=1;
        for(unsigned i=1;i<n;++i)if(x[i]!=x[m-1])x[m++]=x[i];
        return m;
    }
    // Connected components, not just a positive/negative aggregate. A sign
    // region can have two disconnected pieces in a general bilinear map.
    [[nodiscard]] HOST_DEVICE bool components(FoldedPatchRecord& r) const noexcept
    {
        double x[4]{};const unsigned n=cuts(x);r.branch_count=0;
        for(int sign=-1;sign<=1;sign+=2)
        {
            int active=-1;
            for(unsigned i=0;i+1<n;++i)
            {
                const double mid=(x[i]+x[i+1])*0.5;double lo=0,hi=0;
                if(!(mid>x[i]&&mid<x[i+1])) return false;
                if(!section(j,sign,mid,lo,hi)){active=-1;continue;}
                double bl=0,bh=0;
                const bool joined=active>=0 && section(j,sign,x[i],bl,bh);
                if(!joined)
                {
                    if(r.branch_count>=4)return false;
                    active=int(r.branch_count++);
                    auto& b=r.branches[active];b.sign=sign;b.u0=x[i];
                }
                r.branches[active].u1=x[i+1];
            }
        }
        return r.branch_count>=2;
    }
    [[nodiscard]] HOST_DEVICE static int branch_at(const FoldedPatchRecord& r,double u,double v) noexcept
    {
        if(r.flags!=FoldedReady || !(u>=0&&u<=1&&v>=0&&v<=1))return -1;
        const double b=mix(r.jacobian[0],r.jacobian[1],u);
        const double t=mix(r.jacobian[2],r.jacobian[3],u);
        const double jac=mix(b,t,v);
        // Only a sign-uncertainty guard; never clamp a density or move a hit.
        const double bound=(1-u)*(1-v)*::fabs(r.jacobian[0])+u*(1-v)*::fabs(r.jacobian[1])
                          +(1-u)*v*::fabs(r.jacobian[2])+u*v*::fabs(r.jacobian[3]);
        if(::fabs(jac)<=0x1p-46*bound)return -1;
        const int sign=jac>0?1:-1;
        for(unsigned i=0;i<r.branch_count;++i)
            if(r.branches[i].sign==sign && u>=r.branches[i].u0 && u<=r.branches[i].u1)return int(i);
        return -1;
    }
};
} // namespace rainbow
