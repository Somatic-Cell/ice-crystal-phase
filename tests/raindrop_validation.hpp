#pragma once

#include <rainbow/raindrop_paths.hpp>
#include <rainbow/raindrop_settings.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <complex>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace rainbow::tests
{
inline void require(const bool condition,const char* message)
{if(!condition)throw std::runtime_error(message);}
inline double phase_value(const PhaseCycles p){return double(p.turns)+double(p.fraction);}

// 球では独立した quadratic/chord 解析解を使う．陰関数 solver や root bracket を呼ばない．
struct AnalyticSphereIntersector
{
    BoundaryHit intersect(const BoundaryRay& ray) const noexcept
    {
        const auto o=ray.physical_origin.cast<double>(),w=ray.direction.cast<double>();
        const double a=w.dot(w),b=o.dot(w),c=o.dot(o)-1.0;
        const double disc=b*b-a*c;
        BoundaryHit hit{};
        if(disc<0.0)return hit;
        const double root=std::sqrt(disc);
        double t=(-b-root)/a;
        if(t<=double(ray.start_distance)) t=(-b+root)/a;
        if(t<=double(ray.start_distance))return hit;
        const auto p=o.at(w,t);
        hit.code=BoundaryCode::Hit;hit.distance=float(t);hit.position=p.cast<float>();
        hit.outward_normal=p.normalized().cast<float>();
        return hit;
    }
};

// 非球形の比較値．cos(n theta) を直接評価し，dense bracket + bisection を用いる．
// 本番の Chebyshev/interval/Newton は呼ばない．この有限解像度参照は検査 ray の比較専用．
struct CosineReferenceIntersector
{
    const RaindropShape& shape;
    double value(const Vec3T<double> p) const noexcept
    {
        const double r=p.length();
        if(r<1e-14)return -1.0;
        const double mu=std::max(-1.0,std::min(1.0,-p.y/r));
        const double theta=std::acos(mu);
        double radius=1.0;
        for(unsigned n=0;n<8;++n)radius+=double(shape.coefficients[n])*std::cos(double(n)*theta);
        return r-radius;
    }
    BoundaryHit intersect(const BoundaryRay& ray) const noexcept
    {
        BoundaryHit hit{};
        const auto o=ray.physical_origin.cast<double>(),w=ray.direction.cast<double>();
        const double start=ray.start_distance,limit=8.0*shape.outer_radius+4.0;
        double lo=start,fl=value(o.at(w,lo));
        constexpr unsigned samples=4096;
        for(unsigned j=1;j<=samples;++j)
        {
            double hi=start+(limit-start)*double(j)/samples;
            const double fh=value(o.at(w,hi));
            if((fl<0)!=(fh<0))
            {
                for(unsigned iteration=0;iteration<60;++iteration)
                {
                    const double mid=lo+(hi-lo)*0.5,fm=value(o.at(w,mid));
                    if((fm<0)==(fl<0)){lo=mid;fl=fm;}else hi=mid;
                }
                const double t=lo+(hi-lo)*0.5;
                const auto p=o.at(w,t);
                // 法線は陰関数の有限差分で独立に求める．
                constexpr double h=1e-6;
                const Vec3T<double> g{
                    (value(p+Vec3T<double>{h,0,0})-value(p-Vec3T<double>{h,0,0}))/(2*h),
                    (value(p+Vec3T<double>{0,h,0})-value(p-Vec3T<double>{0,h,0}))/(2*h),
                    (value(p+Vec3T<double>{0,0,h})-value(p-Vec3T<double>{0,0,h}))/(2*h)};
                hit.code=BoundaryCode::Hit;hit.position=p.cast<float>();hit.distance=float(t);hit.outward_normal=g.normalized().cast<float>();
                return hit;
            }
            lo=hi;fl=fh;
        }
        return hit;
    }
};

template<class Intersector> inline std::vector<OutgoingVertex> cpu_trace(const RaindropTraceConfig& c,Intersector& intersector)
{
    std::vector<OutgoingVertex> vertices(std::size_t(c.vertex_count())*4);
    for(std::uint32_t i=0;i<c.vertex_count();++i)RaindropPathTracer::trace_vertex(c,i,intersector,vertices.data());
    return vertices;
}

inline void verify_trace(const RaindropTraceConfig& c,const std::vector<OutgoingVertex>& actual,
    const std::vector<OutgoingVertex>& reference)
{
    require(actual.size()==reference.size(),"Trace result size mismatch.");
    std::size_t valid=0,fallback=0;
    double maximum_phase_error=0;
    for(std::size_t i=0;i<actual.size();++i)
    {
        const auto& a=actual[i];const auto& b=reference[i];
        if(a.status!=b.status)
        {
            std::ostringstream s;s<<"Trace status mismatch at "<<i<<": "<<unsigned(a.status)<<" vs "<<unsigned(b.status);throw std::runtime_error(s.str());
        }
        require(a.status==VertexStatus::Valid||a.status==VertexStatus::Miss||a.status==VertexStatus::TotalInternalReflection,
            "Trace contains unresolved/numerical failures.");
        if(a.status!=VertexStatus::Valid)continue;
        ++valid;fallback+=bool(a.diagnostics&IntersectionFallback);
        require(a.position_drop.is_finite()&&a.direction_drop.is_finite()&&a.basis_x.is_finite(),"Nonfinite vertex.");
        require(std::abs(a.direction_drop.length()-1)<3e-6f,"Nonunit direction.");
        require(std::abs(a.direction_drop.dot(a.basis_x))<4e-6f,"Polarization frame is not transverse.");
        require(std::abs(a.basis_x.length()-1)<4e-6f,"Nonunit polarization basis.");
        require(std::abs(c.shape.implicit_value(a.position_drop))<8e-6f,"Hit is off the surface.");
        require((a.position_drop-b.position_drop).length()<3e-4f,"Exit position discrepancy.");
        require((a.direction_drop-b.direction_drop).length()<3e-4f,"Exit direction discrepancy.");
        const double phase_error=std::abs(phase_value(a.optical_cycles)-phase_value(b.optical_cycles));
        maximum_phase_error=std::max(maximum_phase_error,phase_error);
        // これは vertex の幾何誤差を含む検査閾値．最終 fringe/LUT 誤差の保証ではない．
        require(a.optical_cycles.is_valid()&&phase_error<0.06,"Optical-path cycles discrepancy.");
        const double na=a.field.squared_norm(),nb=b.field.squared_norm();
        require(std::isfinite(na)&&std::abs(na-nb)<2e-4*(1+std::abs(nb)),"Field amplitude discrepancy.");
    }
    require(valid>0,"No valid outgoing vertices.");
    std::cout<<"trace: "<<valid<<" valid, "<<fallback<<" intersection-fallback vertices, max cycle error "<<maximum_phase_error<<'\n';
}

inline void verify_sphere_optical_lengths(const RaindropTraceConfig& c,const std::vector<OutgoingVertex>& vertices)
{
    // 独立な式: R: d+d'-2cos(i); T R^(p-1) T: d+d'-2cos(i)+2 p n cos(r)．
    const double scale=double(c.radius_mm)*1e6/double(c.wavelength_nm);
    for(std::uint32_t i=0;i<c.vertex_count();++i)
    {
        const float uf=(2.0f*float(i%c.grid_width)/float(c.grid_width-1)-1.0f)*c.grid_half_extent;
        const float vf=(2.0f*float(i/c.grid_width)/float(c.grid_height-1)-1.0f)*c.grid_half_extent;
        const double b2=double(uf)*uf+double(vf)*vf;
        if(b2>=0.98)continue; // grazing を解析値精度検査とは分離（上の比較検査では除外しない）．
        const double ci=std::sqrt(1-b2),ct=std::sqrt(1-b2/(double(c.interior_index)*c.interior_index));
        for(unsigned p=0;p<4;++p)
        {
            const auto& v=vertices[std::size_t(p)*c.vertex_count()+i];
            require(v.status==VertexStatus::Valid,"Sphere should emit every path family.");
            const double expected=(double(c.reference_distance)+c.outgoing_reference_distance-2*ci
                +2*double(p)*double(c.interior_index)*ct)*scale;
            require(std::abs(phase_value(v.optical_cycles)-expected)<0.06,"Sphere analytic optical length mismatch.");
            const double incident_angle=std::asin(std::sqrt(b2));
            const double refracted_angle=std::asin(std::sqrt(b2)/double(c.interior_index));
            constexpr double pi=3.14159265358979323846;
            const double deviation=p==0?pi-2*incident_angle:
                double(p-1)*pi+2*incident_angle-2*double(p)*refracted_angle;
            require(std::abs(double(c.incident_direction.dot(v.direction_drop))-std::cos(deviation))<5e-5,
                "Sphere analytic deflection angle mismatch.");
        }
    }
}
inline void verify_fresnel()
{
    const Vec3 normal{0,1,0};
    for(const float cosine:{1.0f,0.95f,0.8f,0.5f,0.1f})
    {
        const Vec3 w{std::sqrt(1-cosine*cosine),-cosine,0};
        const auto f=TransverseFrame::from_direction(w);
        for(unsigned pol=0;pol<2;++pol)
        {
            PolarizedRay ray{w,f,pol?Field32{{0,0},{1,0}}:Field32{{1,0},{0,0}}};
            const auto s=DielectricInterface::split(ray,normal,1.0f,1.3314f);
            require(s.is_valid&&s.has_transmission,"Air/water Fresnel split failed.");
            const double ct=-s.transmitted.direction.dot(normal);
            const double power=s.reflected.field.squared_norm()+1.3314*ct/cosine*s.transmitted.field.squared_norm();
            require(std::abs(power-1)<2e-6,"Interface power conservation failure.");
            // ベクトルとして接線 E_i + E_r = E_t を検査する．power だけでは符号誤りを見逃す．
            const auto physical_field=[](const PolarizedRay& p)
            {return p.frame.e0*p.field.x.real+p.frame.e1*p.field.y.real;};
            const Vec3 electric_difference=physical_field(ray)+physical_field(s.reflected)-physical_field(s.transmitted);
            const Vec3 tangential=electric_difference-normal*electric_difference.dot(normal);
            require(tangential.length()<3e-6f,"Tangential electric-field continuity/sign mismatch.");
        }
    }
    // TIR: 実振幅への丸めではなく，非零の虚部と |r|=1 を確認する．
    PolarizedRay ray{{0.8f,-0.6f,0},{{0,0,1},{-0.6f,-0.8f,0}},{{1,0},{1,0}}};
    const auto tir=DielectricInterface::split(ray,normal,1.5f,1.0f);
    require(tir.is_valid&&!tir.has_transmission,"TIR branch mismatch.");
    require(std::abs(tir.reflected.field.x.squared_norm()-1)<1e-6f&&std::abs(tir.reflected.field.y.squared_norm()-1)<1e-6f,"TIR coefficient modulus mismatch.");
    require(std::abs(tir.reflected.field.x.imag)>0.1f,"TIR complex phase was lost.");
    const double ci=-double(ray.direction.y);
    const double beta=std::sqrt(1.5*1.5*double(ray.direction.x)*double(ray.direction.x)-1.0);
    const auto rs=std::polar(1.0,-2.0*std::atan2(beta,1.5*ci));
    const auto rp=std::polar(1.0,-2.0*std::atan2(1.5*beta,ci));
    require(std::abs(double(tir.reflected.field.x.real)-rs.real())<2e-6
        &&std::abs(double(tir.reflected.field.x.imag)-rs.imag())<2e-6
        &&std::abs(double(tir.reflected.field.y.real)-rp.real())<2e-6
        &&std::abs(double(tir.reflected.field.y.imag)-rp.imag())<2e-6,"TIR phase sign/value mismatch.");
}
inline void verify_shape_table()
{
    for(float radius:{0.005f,0.4f,0.7f,1.0f,1.5f,2.0f,2.5f,3.0f})
    {
        RaindropShape shape{};require(RaindropShape::try_make(radius,false,shape),"Shape setup failed.");
        CosineReferenceIntersector reference{shape};
        for(unsigned i=0;i<=512;++i)
        {
            const double theta=3.14159265358979323846*double(i)/512;
            double slope;
            const double r=shape.radius_and_derivative(std::cos(theta),slope);
            const Vec3T<double> p{r*std::sin(theta),-r*std::cos(theta),0};
            require(std::abs(reference.value(p))<2e-14,"Chebyshev/cosine shape mismatch.");
            require(r>=shape.inner_radius&&r<=shape.outer_radius,"Nonconservative shape radius bound.");
        }
    }
}
inline void verify_vertical_orientation()
{
    // +y は上．r(theta) の theta=0 を下向きにして初めて下側が扁平になる．
    // 極の曲率 k=(r-r_theta_theta)/r^2 を独立に Table I の cosine 係数から計算する．
    RaindropShape shape{};
    require(RaindropShape::try_make(3.0f,false,shape),"Orientation shape setup failed.");
    double bottom_radius=1.0,top_radius=1.0,bottom_second=0,top_second=0;
    for(unsigned n=0;n<8;++n)
    {
        const double c=shape.coefficients[n],sign=(n%2)?-1.0:1.0;
        bottom_radius+=c;top_radius+=sign*c;
        bottom_second-=double(n*n)*c;top_second-=double(n*n)*sign*c;
    }
    const double bottom_curvature=(bottom_radius-bottom_second)/(bottom_radius*bottom_radius);
    const double top_curvature=(top_radius-top_second)/(top_radius*top_radius);
    require(bottom_curvature<top_curvature*0.25,"Drop should have a flatter bottom, not a flatter top.");
    require(std::abs(shape.implicit_value(Vec3T<double>{0,-bottom_radius,0}))<1e-14,
        "Bottom pole does not match the downward polar-angle convention.");
    require(std::abs(shape.implicit_value(Vec3T<double>{0,top_radius,0}))<1e-14,
        "Top pole does not match the downward polar-angle convention.");
}
inline void verify_compensated_optical_path()
{
    OpticalPathAccumulator path=OpticalPathAccumulator::make(3.0f,380.0f);
    double expected=0;
    for(unsigned i=0;i<200;++i)
    {
        const float length=0.013f+float(i%7)*0.071f;
        require(path.add_distance({length,0},1.3314f),"Optical path update failed.");
        expected+=double(length)*double(1.3314f)*(double(3.0f)*1e6/double(380.0f));
    }
    require(std::abs(phase_value(path.phase)-expected)<2e-5,"Compensated optical-path arithmetic mismatch.");
}
}
