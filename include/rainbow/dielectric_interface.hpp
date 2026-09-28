#pragma once

#include <rainbow/field32.hpp>
#include <rainbow/raindrop_intersection.hpp>

namespace rainbow
{

struct PolarizedRay
{
    Vec3 direction;
    TransverseFrame frame;
    Field32 field;
};
struct InterfaceBranches
{
    PolarizedRay reflected;
    PolarizedRay transmitted;
    bool has_transmission=false;
    bool is_valid=false;
    std::uint32_t diagnostics=0;
};

// 電場振幅係数を使う．反射・透過の選択，ロシアンルーレット，振幅による枝刈りはしない．
// 固定規約: s = normalize(w_i x n), p_i = w_i x s,
// p_r = w_r x s, p_t = w_t x s．n は入射媒質側を向く法線．
struct DielectricInterface
{
    HOST_DEVICE static InterfaceBranches split(
        const PolarizedRay& incoming,const Vec3 oriented_normal,
        const float n_i,const float n_t) noexcept
    {
        InterfaceBranches output{};
        const Vec3 wi=incoming.direction;
        const Vec3 n=oriented_normal;
        float cosine=-wi.dot(n);
        if(!(cosine>=-0x1p-19f && cosine<=1.0f+0x1p-19f && n_i>0 && n_t>0)) return output;
        cosine=::fminf(1.0f,::fmaxf(0.0f,cosine));
        Vec3 s=wi.cross(n);
        const float sin2=s.dot(s);
        if(sin2>0x1p-40f) s=s/::sqrtf(sin2);
        else s=(incoming.frame.e0-wi*wi.dot(incoming.frame.e0)).normalized();
        const Vec3 pi=wi.cross(s).normalized();
        const Complex32 es=incoming.field.x*incoming.frame.e0.dot(s)
            +incoming.field.y*incoming.frame.e1.dot(s);
        const Complex32 ep=incoming.field.x*incoming.frame.e0.dot(pi)
            +incoming.field.y*incoming.frame.e1.dot(pi);

        const float eta=n_i/n_t;
        float discriminant=::fmaf(-(eta*eta),sin2,1.0f);
        // 分岐が変わり得る critical-angle 近傍だけ double で再評価する．
        if(::fabsf(discriminant)<0x1p-15f)
        {
            const auto wd=wi.cast<double>().normalized(),nd=n.cast<double>().normalized();
            const double ratio=double(n_i)/double(n_t);
            const auto sd=wd.cross(nd);
            discriminant=static_cast<float>(::fma(-ratio*ratio,sd.dot(sd),1.0));
            cosine=static_cast<float>(-wd.dot(nd));
            output.diagnostics|=SnellFallback;
        }
        const Vec3 wr=wi.at(n,2.0f*cosine).normalized();
        output.reflected.direction=wr;
        output.reflected.frame={s,wr.cross(s).normalized()};
        Complex32 rs{},rp{};
        if(discriminant<0.0f)
        {
            const float beta=::sqrtf(-discriminant);
            rs=tir_coefficient(n_i*cosine,n_t*beta);
            rp=tir_coefficient(n_t*cosine,n_i*beta);
        }
        else
        {
            const float ct=::sqrtf(discriminant);
            const float ds=n_i*cosine+n_t*ct;
            const float dp=n_t*cosine+n_i*ct;
            if(!(ds>0.0f && dp>0.0f)) return output;
            rs={(n_i*cosine-n_t*ct)/ds,0};
            rp={(n_t*cosine-n_i*ct)/dp,0};
            const float ts=2.0f*n_i*cosine/ds;
            const float tp=2.0f*n_i*cosine/dp;
            const Vec3 wt=(wi*eta+n*(eta*cosine-ct)).normalized();
            output.transmitted={wt,{s,wt.cross(s).normalized()},{es*ts,ep*tp}};
            output.has_transmission=true;
        }
        output.reflected.field={es*rs,ep*rp};
        output.is_valid=wr.is_finite() && (!output.has_transmission||output.transmitted.direction.is_finite());
        return output;
    }
private:
    HOST_DEVICE static Complex32 tir_coefficient(const float a,const float b) noexcept
    {
        // (a-i b)/(a+i b)．正の符号規約 exp(i kz-i wt) に対する全反射位相．
        const float scale=::fmaxf(::fabsf(a),::fabsf(b));
        const float x=a/scale,y=b/scale,den=::fmaf(x,x,y*y);
        return {::fmaf(x,x,-y*y)/den,(-2.0f*x*y)/den};
    }
};
}
