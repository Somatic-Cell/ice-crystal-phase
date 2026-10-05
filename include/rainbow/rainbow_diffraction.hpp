#pragma once
#include <rainbow/wave_optics_data.hpp>
#include <cmath>

namespace rainbow
{
// Explicit discrete realization of Sec.4.1.2, not the unavailable author code.
// Table-II bandwidths; primary/secondary support transition labels; spherical
// Gaussian on TOTAL coherent intensity; solid-angle quadrature; 4-sigma support;
// activation is tapered from 3 to 4 sigma. Unavailable samples are NOT filled.
struct RainbowDiffraction
{
    // Closed-cap comparisons tolerate only FP64 roundoff at the cutoff.
    static constexpr double angular_boundary_tolerance=0x1p-44;
    static constexpr double pi=3.141592653589793238462643383279502884;
    [[nodiscard]] static bool table_sigma_degrees(float radius,double& output) noexcept
    {
        constexpr float r[10]={0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f};
        constexpr double s[10]={0.70,0.45,0.30,0.25,0.22,0.20,0.18,0.17,0.16,0.15};
        if(!(radius>=r[0] && radius<=r[9])) return false;
        for(unsigned j=0;j<10;++j)
        {
            if(radius==r[j]){output=s[j];return true;}
            if(j && radius<r[j])
            {output=::fma((double(radius)-r[j-1])/(double(r[j])-r[j-1]),s[j]-s[j-1],s[j-1]);return true;}
        }
        return false;
    }
    [[nodiscard]] HOST_DEVICE static bool valid_config(const DiffractionParams& p) noexcept
    {
        const auto& c=p.config;
        return p.theta_count>=2u && p.phi_count>=1u
            && p.theta_count<=32767u && p.phi_count<=32767u
            && c.primary_sigma_rad>0 && finite(c.primary_sigma_rad)
            && c.support_sigma>=2 && c.support_sigma<=8
            && c.inner_blend_sigma>=0 && c.inner_blend_sigma<c.support_sigma
            && 2*c.primary_sigma_rad*c.support_sigma<pi
            && c.transition_contrast>=0 && c.transition_contrast<1;
    }
    [[nodiscard]] HOST_DEVICE static RainbowTransition detect(const DiffractionParams& p,std::uint32_t i) noexcept
    {
        RainbowTransition r{};
        if(!valid_config(p) || !p.focal || i>=p.theta_count*p.phi_count)
        {r.flags=TransitionUnknown;return r;}
        if(i/p.phi_count+1u==p.theta_count)return r;
        const auto& a=p.focal[i];const auto& b=p.focal[i+p.phi_count];
        if(!sample_valid(a)||!sample_valid(b)){r.flags=TransitionUnknown;return r;}
        const double x=a.intensity_s+a.intensity_p,y=b.intensity_s+b.intensity_p;
        const double peak=::fmax(x,y);
        if(!(peak>0) || !finite(peak) || ::fabs(y-x)/peak<p.config.transition_contrast)return r;
        // Primary: increasing theta first enters TRT support. Secondary: leaves TRRT.
        // Gate by family so ordinary fringes and the TT cutoff are NOT all blurred.
        if(a.family_hits[2]==0 && b.family_hits[2]>0 && b.family_incoherent[2]>0)r.flags|=TransitionPrimary;
        if(a.family_hits[3]>0 && b.family_hits[3]==0 && a.family_incoherent[3]>0)r.flags|=TransitionSecondary;
        return r;
    }
    [[nodiscard]] HOST_DEVICE static DiffractionResult filter(const DiffractionParams& p,std::uint32_t index) noexcept
    {
        DiffractionResult out{};
        if(!valid_config(p)||!p.focal||!p.transitions||index>=p.theta_count*p.phi_count)
        {out.flags=DiffractionInvalidData;invalidate(out);return out;}
        const auto& center=p.focal[index];
        if(!sample_valid(center)){out.flags=DiffractionInputUnavailable;invalidate(out);return out;}
        out.intensity_s=center.intensity_s;out.intensity_p=center.intensity_p;
        const int nr=int(p.theta_count),nc=int(p.phi_count),row=int(index/p.phi_count),col=int(index%p.phi_count);
        const double dt=pi/double(nr),dp=2*pi/double(nc),t=(row+0.5)*dt;
        const double support=p.config.support_sigma;
        const double search_radius=2*p.config.primary_sigma_rad*support;
        const int lo=maxi(0,int(::floor((t-search_radius)/dt))-1);
        const int hi=mini(nr-2,int(::ceil((t+search_radius)/dt)));
        double best=support+1;
        // Find nearest normalized-distance transition; overlap uses that label's
        // bandwidth (deterministic tie order). This selection is a convention.
        for(int rr=lo;rr<=hi;++rr)
        {
            const double tr=(rr+1.0)*dt;
            const int width=phi_radius(t,tr,search_radius,dp,nc);
            if(width<0)continue;
            const int count=mini(nc,2*width+1),start=count==nc?0:col-width;
            for(int j=0;j<count;++j)
            {
                const int cc=wrap(start+j,nc);
                const double d=angular_distance(t,tr,(cc-col)*dp);
                if(d>search_radius+angular_boundary_tolerance)continue;
                const auto flag=p.transitions[std::uint32_t(rr)*p.phi_count+std::uint32_t(cc)].flags;
                if(flag&TransitionUnknown)out.flags|=DiffractionDetectionUnavailable;
                for(unsigned kind=TransitionPrimary;kind<=TransitionSecondary;kind*=2u)
                {
                    if(!(flag&kind))continue;
                    const double sigma=p.config.primary_sigma_rad*(kind==TransitionSecondary?2:1);
                    const double distance=d/sigma;
                    if(distance<support && (distance<best-1e-12
                        || (::fabs(distance-best)<=1e-12 && kind<out.transition_kind)))
                    {best=distance;out.sigma_rad=sigma;out.transition_kind=kind;}
                }
            }
        }
        if(out.flags&DiffractionDetectionUnavailable){invalidate(out);return out;}
        if(out.transition_kind==0)return out; // unchanged, no artificial global blur
        const double sigma=out.sigma_rad;
        if(dt>0.5*sigma || ::sin(t)*dp>0.5*sigma)out.flags|=DiffractionUnderresolved;
        out.blend=1;
        if(best>p.config.inner_blend_sigma)
        {
            const double x=(best-p.config.inner_blend_sigma)/(support-p.config.inner_blend_sigma);
            out.blend=1-x*x*(3-2*x);
        }
        const double radius=support*sigma;
        const int lower=maxi(0,int(::floor((t-radius)/dt))-1);
        const int upper=mini(nr-1,int(::ceil((t+radius)/dt)));
        Sum sum_s{},sum_p{},weights{};
        for(int rr=lower;rr<=upper;++rr)
        {
            const double tr=(rr+0.5)*dt;
            const int width=phi_radius(t,tr,radius,dp,nc);
            if(width<0)continue;
            const int count=mini(nc,2*width+1),start=count==nc?0:col-width;
            const double area=2*::sin(tr)*::sin(0.5*dt)*dp;
            for(int j=0;j<count;++j)
            {
                const int cc=wrap(start+j,nc);
                const double d=angular_distance(t,tr,(cc-col)*dp);
                if(d>radius+angular_boundary_tolerance)continue;
                const auto& sample=p.focal[std::uint32_t(rr)*p.phi_count+std::uint32_t(cc)];
                if(!sample_valid(sample)){out.flags|=DiffractionStencilUnavailable;continue;}
                const double weight=::exp(-0.5*(d/sigma)*(d/sigma))*area;
                weights.add(weight);sum_s.add(weight*sample.intensity_s);sum_p.add(weight*sample.intensity_p);
            }
        }
        if(out.flags&DiffractionStencilUnavailable){invalidate(out);return out;}
        const double w=weights.value();
        if(!(w>0) || !finite(w)){out.flags|=DiffractionInvalidData;invalidate(out);return out;}
        // Constant-preserving local normalization, NOT global energy renormalization.
        out.intensity_s=::fma(out.blend,sum_s.value()/w-center.intensity_s,center.intensity_s);
        out.intensity_p=::fma(out.blend,sum_p.value()/w-center.intensity_p,center.intensity_p);
        if(!(out.intensity_s>=0 && out.intensity_p>=0) || !finite(out.intensity_s+out.intensity_p))
        {out.flags|=DiffractionInvalidData;invalidate(out);return out;}
        out.flags|=DiffractionFiltered;return out;
    }
    [[nodiscard]] HOST_DEVICE static double angular_distance(double a,double b,double phi) noexcept
    {
        const double u=::sin(0.5*(a-b)),v=::sin(0.5*phi);
        const double h=::fma(::sin(a)*::sin(b),v*v,u*u);
        return 2*::asin(::sqrt(::fmin(1.0,::fmax(0.0,h))));
    }
private:
    struct Sum
    {
        double s=0,c=0;
        HOST_DEVICE void add(double x) noexcept
        {const double n=s+x;c+=::fabs(s)>=::fabs(x)?(s-n)+x:(x-n)+s;s=n;}
        [[nodiscard]] HOST_DEVICE double value()const noexcept{return s+c;}
    };
    [[nodiscard]] HOST_DEVICE static bool finite(double x) noexcept
    {return x>=-0x1.fffffffffffffp1023 && x<=0x1.fffffffffffffp1023;}
    [[nodiscard]] HOST_DEVICE static bool sample_valid(const FocalOpticalResult& r) noexcept
    {return r.valid() && r.intensity_s>=0 && r.intensity_p>=0 && finite(r.intensity_s+r.intensity_p);}
    HOST_DEVICE static void invalidate(DiffractionResult& r)noexcept
    {r.intensity_s=r.intensity_p=patch_optical_nan;}
    HOST_DEVICE static int maxi(int a,int b)noexcept{return a>b?a:b;}
    HOST_DEVICE static int mini(int a,int b)noexcept{return a<b?a:b;}
    HOST_DEVICE static int wrap(int a,int n)noexcept{const int r=a%n;return r<0?r+n:r;}
    // Conservative longitude window for a spherical cap. Exact angular test
    // follows; columns are never duplicated even when the cap includes a pole.
    HOST_DEVICE static int phi_radius(double t,double r,double radius,double dp,int nc)noexcept
    {
        if(::fabs(t-r)>radius+angular_boundary_tolerance)return -1;
        const double a=::sin(t)*::sin(r),b=::cos(t)*::cos(r);
        if(!(a>0))return nc;
        const double x=(::cos(radius)-b)/a;
        if(x<=-1)return nc;
        if(x>=1)return 1; // boundary roundoff: exact distance test rejects extras
        return mini(nc,int(::ceil(::acos(x)/dp))+1);
    }
};
} // namespace rainbow
