#pragma once

#include <rainbow/vec3.hpp>
#include <cmath>
#include <cstdint>

namespace rainbow
{

// Sadeghi 2012 Eq.(1), Table I (p.5)．座標は等体積球半径 a で割った無次元座標．
// +y を上，gravity を -y とし，表の極角 theta は -y 軸から測る．
// この規約で負の奇数次係数が下側の扁平化を表す．c0 も和に含める．
struct RaindropShape
{
    float coefficients[8] = {};
    float inner_radius = 1.0f;
    float outer_radius = 1.0f;

    // 半径は diameter ではない．表の範囲外へ勝手に外挿しない．
    // a<=0.4mm は表の球の行．sphere=true は比較検証のため明示的に球を指定する．
    static bool try_make(const float radius_mm, const bool sphere, RaindropShape& result) noexcept
    {
        if(!(radius_mm>0.0f && radius_mm<=3.0f)) return false;
        constexpr float radii[6] = {0.4f,1.0f,1.5f,2.0f,2.5f,3.0f};
        constexpr float rows[6][8] = {
            {0,0,0,0,0,0,0,0},
            {-0.0131f,-0.0120f,-0.0376f,-0.0096f,-0.0004f,0.0015f,0.0005f,0},
            {-0.0282f,-0.0230f,-0.0779f,-0.0175f,0.0021f,0.0046f,0.0011f,-0.0006f},
            {-0.0458f,-0.0335f,-0.1211f,-0.0227f,0.0083f,0.0089f,0.0012f,-0.0021f},
            {-0.0644f,-0.0416f,-0.1629f,-0.0246f,0.0176f,0.0131f,0.0002f,-0.0044f},
            {-0.0840f,-0.0480f,-0.2034f,-0.0237f,0.0297f,0.0166f,-0.0021f,-0.0072f}
        };
        RaindropShape shape{};
        if(!sphere && radius_mm>radii[0])
        {
            unsigned upper=1;
            while(upper<5 && radius_mm>radii[upper]) ++upper;
            const float t=(radius_mm-radii[upper-1])/(radii[upper]-radii[upper-1]);
            for(unsigned k=0;k<8;++k)
                shape.coefficients[k]=::fmaf(t,rows[upper][k]-rows[upper-1][k],rows[upper-1][k]);
        }
        // |T_n(mu)|<=1．幾何そのものの近似でなく，保守的な包囲球だけに使う．
        float absolute_sum=0.0f;
        for(unsigned k=1;k<8;++k)
            absolute_sum=::nextafterf(absolute_sum+::fabsf(shape.coefficients[k]),INFINITY);
        const float center=1.0f+shape.coefficients[0];
        shape.inner_radius=::nextafterf(center-absolute_sum-0x1p-20f,-INFINITY);
        shape.outer_radius=::nextafterf(center+absolute_sum+0x1p-20f,INFINITY);
        result=shape;
        return true;
    }

    // cos(n theta)=T_n(cos theta)．acos や極での sin(theta) による除算は不要．
    template<class T> HOST_DEVICE T radius_and_derivative(const T mu, T& derivative) const noexcept
    {
        T previous=T(1), current=mu, d_previous=T(0), d_current=T(1);
        T radius=T(1)+T(coefficients[0])+T(coefficients[1])*mu;
        derivative=T(coefficients[1]);
        for(unsigned n=2;n<8;++n)
        {
            const T next=T(2)*mu*current-previous;
            const T d_next=T(2)*current+T(2)*mu*d_current-d_previous;
            radius+=T(coefficients[n])*next;
            derivative+=T(coefficients[n])*d_next;
            previous=current; current=next; d_previous=d_current; d_current=d_next;
        }
        return radius;
    }
    template<class T> HOST_DEVICE T implicit_value(const Vec3T<T> p) const noexcept
    {
        const T r=p.length();
        if(r<T(inner_radius)*T(0.25)) return r-T(inner_radius); // 内部の符号だけを返す領域．root は存在しない．
        T derivative;
        return r-radius_and_derivative(-p.y/r,derivative);
    }
    template<class T> HOST_DEVICE Vec3T<T> gradient(const Vec3T<T> p) const noexcept
    {
        const T r=p.length();
        const Vec3T<T> u=p/r;
        T derivative;
        const T mu=-u.y;
        static_cast<void>(radius_and_derivative(mu,derivative));
        return u-(Vec3T<T>{T(0),T(-1),T(0)}-u*mu)*(derivative/r);
    }
};
static_assert(sizeof(RaindropShape)==40);
}
