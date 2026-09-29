#pragma once

#include <rainbow/raindrop_trace_data.hpp>
#include <cstdint>
#include <cmath>
#include <limits>
#include <type_traits>

namespace rainbow
{

enum class BoundaryCode : std::uint32_t { Miss=0, Hit=1, Unresolved=2 };

struct BoundaryRay
{
    Vec3 physical_origin;
    Vec3 direction;
    float start_distance = 0.0f; // 数値的 self-hit 除外．光路に加算する距離ではない．
    bool starts_inside = false;
};
struct BoundaryHit
{
    BoundaryCode code = BoundaryCode::Miss;
    Vec3 position;
    Vec3 outward_normal;
    float distance = 0.0f; // physical_origin から測る．start_distance を引かない．
    std::uint32_t diagnostics = 0;
};

namespace detail
{
// 型ごとの正の無限大．namespace scope の constexpr 値として用意する．
// GPU 側ではこの値を使用し，std::numeric_limits の関数を直接呼ばない．
template<class T>
inline constexpr T positive_infinity_v =
    std::numeric_limits<T>::infinity();

template<class T>
HOST_DEVICE T next_down(const T x) noexcept
{
    static_assert(
        std::is_same_v<T, float> || std::is_same_v<T, double>);

    if constexpr(std::is_same_v<T, float>)
    {
        return ::nextafterf(x, -positive_infinity_v<float>);
    }
    else
    {
        return ::nextafter(x, -positive_infinity_v<double>);
    }
}

template<class T>
HOST_DEVICE T next_up(const T x) noexcept
{
    static_assert(
        std::is_same_v<T, float> || std::is_same_v<T, double>);

    if constexpr(std::is_same_v<T, float>)
    {
        return ::nextafterf(x, positive_infinity_v<float>);
    }
    else
    {
        return ::nextafter(x, positive_infinity_v<double>);
    }
}
template<class T> HOST_DEVICE T root_sqrt(const T x) noexcept
{
    if constexpr(std::is_same_v<T,float>) return ::sqrtf(x);
    else return ::sqrt(x);
}
template<class T> HOST_DEVICE T abs_value(const T x) noexcept { return x<T(0)?-x:x; }
template<class T> HOST_DEVICE T max_value(const T a,const T b) noexcept { return a>b?a:b; }
template<class T> HOST_DEVICE T min_value(const T a,const T b) noexcept { return a<b?a:b; }

// root 候補の除外専用の外向き丸め interval．係数は格納された FP32 値を正確な入力とする．
// fast-math / FTZ / 式の再結合を有効にしないこと．0 が含まれる区間を miss と決めつけない．
template<class T> struct RootInterval
{
    T lo,hi;
    HOST_DEVICE static RootInterval point(const T x) noexcept { return {x,x}; }
    HOST_DEVICE RootInterval operator+(const RootInterval b) const noexcept
    { return {next_down(lo+b.lo),next_up(hi+b.hi)}; }
    HOST_DEVICE RootInterval operator-() const noexcept { return {-hi,-lo}; }
    HOST_DEVICE RootInterval operator-(const RootInterval b) const noexcept { return *this+(-b); }
    HOST_DEVICE RootInterval operator*(const RootInterval b) const noexcept
    {
        const T a=lo*b.lo,c=lo*b.hi,d=hi*b.lo,e=hi*b.hi;
        return {next_down(min_value(min_value(a,c),min_value(d,e))),
                next_up(max_value(max_value(a,c),max_value(d,e)))};
    }
    HOST_DEVICE RootInterval operator/(const RootInterval b) const noexcept
    {
        if(b.lo<=T(0) && b.hi>=T(0)) return {-positive_infinity_v<T>,positive_infinity_v<T>};
        return *this * RootInterval{next_down(T(1)/b.hi),next_up(T(1)/b.lo)};
    }
    HOST_DEVICE RootInterval squared() const noexcept
    {
        const T a=lo*lo,b=hi*hi;
        return {lo<=T(0)&&hi>=T(0)?T(0):max_value(T(0),next_down(min_value(a,b))),
                next_up(max_value(a,b))};
    }
    HOST_DEVICE RootInterval square_root() const noexcept
    { return {max_value(T(0),next_down(root_sqrt(max_value(T(0),lo)))),next_up(root_sqrt(max_value(T(0),hi)))}; }
    HOST_DEVICE bool excludes_zero() const noexcept { return lo>T(0) || hi<T(0); }
};

template<class T> struct RootBounds
{
    RootInterval<T> value;
    RootInterval<T> derivative;
};

template<class T> HOST_DEVICE RootBounds<T> bound_surface(
    const RaindropShape& shape,const Vec3T<T> origin,const Vec3T<T> direction,
    const T lower,const T upper) noexcept
{
    using I=RootInterval<T>;
    const I t{lower,upper};
    const I px=I::point(origin.x)+t*I::point(direction.x);
    const I py=I::point(origin.y)+t*I::point(direction.y);
    const I pz=I::point(origin.z)+t*I::point(direction.z);
    const I r=(px.squared()+py.squared()+pz.squared()).square_root();
    if(r.hi<T(shape.inner_radius)) return {{-T(1),-T(1)},{-positive_infinity_v<T>,positive_infinity_v<T>}};
    if(r.lo>T(shape.outer_radius)) return {{T(1),T(1)},{-positive_infinity_v<T>,positive_infinity_v<T>}};
    if(r.lo<=T(0)) return {{-positive_infinity_v<T>,positive_infinity_v<T>},{-positive_infinity_v<T>,positive_infinity_v<T>}};
    I mu=(-py)/r;
    mu.lo=max_value(mu.lo,-T(1)); mu.hi=min_value(mu.hi,T(1));
    I previous=I::point(T(1)),current=mu,dp=I::point(T(0)),dc=I::point(T(1));
    I radius=I::point(T(1))+I::point(T(shape.coefficients[0]))
        +I::point(T(shape.coefficients[1]))*mu;
    I slope=I::point(T(shape.coefficients[1]));
    const I two=I::point(T(2));
    for(unsigned n=2;n<8;++n)
    {
        const I next=two*mu*current-previous;
        const I dn=two*current+two*mu*dc-dp;
        // ゼロ係数も省略可能だが，shape 自体を近似する処理は入れない．
        if(shape.coefficients[n]!=0.0f)
        { radius=radius+I::point(T(shape.coefficients[n]))*next; slope=slope+I::point(T(shape.coefficients[n]))*dn; }
        previous=current; current=next; dp=dc; dc=dn;
    }
    const I dr=(px*I::point(direction.x)+py*I::point(direction.y)+pz*I::point(direction.z))/r;
    const I dm=(I::point(-direction.y)-mu*dr)/r;
    return {r-radius,dr-slope*dm};
}

// 終了理由を数値計算の結果と分離する．診断は判定・許容誤差を変更しない．
enum class RootStopReason : std::uint32_t
{
    None = 0,
    InsideStartNotInside = 1,
    VisitLimit = 2,
    BracketIterationLimit = 3,
    WidthLimit = 4,
    MidpointStagnation = 5,
    StackLimit = 6
};

template<class T> struct RootSolverDiagnostics
{
    BoundaryCode result = BoundaryCode::Miss;
    RootStopReason reason = RootStopReason::None;
    std::uint32_t visits = 0;
    std::uint32_t iterations = 0;
    std::uint32_t pending_segments = 0;
    T lower = T(0), upper = T(0), tolerance = T(0);
    T f_lower = T(0), f_upper = T(0);
    RootBounds<T> bounds{};
    bool has_bounds = false;
    bool has_endpoint_values = false;

    HOST_DEVICE void stop(
        const RootStopReason why, const unsigned visit_count,
        const unsigned iteration_count, const unsigned pending_count,
        const T lo, const T hi, const T tol) noexcept
    {
        reason = why;
        visits = visit_count;
        iterations = iteration_count;
        pending_segments = pending_count;
        lower = lo;
        upper = hi;
        tolerance = tol;
    }
};

// 任意の診断出力．GPU 出力 ABI (OutgoingVertex) には含めない．
struct BoundarySolveDiagnostics
{
    RootSolverDiagnostics<float> fp32{};
    RootSolverDiagnostics<double> fp64{};
    bool used_fp64 = false;
};

template<class T> struct SurfaceRoot { BoundaryCode code; T distance; };

// 陰関数そのものへの bracket + safeguarded Newton．三角形化・楕円体近似はしない．
// 左区間を先に処理し，単調性を interval derivative で確認した bracket だけを解く．
// root を除外できないまま解像度/計算上限に達したら Unresolved を返す（miss へ変換しない）．
template<class T> HOST_DEVICE SurfaceRoot<T> solve_surface(
    const RaindropShape& shape,const BoundaryRay& ray,
    RootSolverDiagnostics<T>* diagnostic = nullptr) noexcept
{
    if(diagnostic != nullptr) *diagnostic = {};
    const Vec3T<T> o=ray.physical_origin.template cast<T>();
    const Vec3T<T> w=ray.direction.template cast<T>();
    const T start=T(ray.start_distance);
    const T limit=T(8)*T(shape.outer_radius)+T(4); // 入射面 d=2R, offset 付き内部点を十分包囲．
    const T tolerance=std::is_same_v<T,float>?T(0x1p-22):T(0x1p-44);

    if(ray.starts_inside)
    {
        const T f=shape.implicit_value(o.at(w,start));
        // ごく浅い grazing chord を self-hit epsilon で黙って飛ばさない．
        if(!(f<T(0)))
        {
            if(diagnostic != nullptr)
            {
                diagnostic->stop(RootStopReason::InsideStartNotInside,0,0,0,start,start,tolerance);
                diagnostic->f_lower = diagnostic->f_upper = f;
                diagnostic->has_endpoint_values = true;
            }
            return {BoundaryCode::Unresolved,T(0)};
        }
    }

    struct Segment { T lo,hi; };
    Segment stack[64];
    unsigned size=0,visits=0;
    stack[size++]={start,limit};
    constexpr unsigned visit_limit=16384;
    while(size!=0)
    {
        if(++visits>visit_limit)
        {
            if(diagnostic != nullptr)
            {
                const Segment pending = stack[size-1];
                const T middle = pending.lo+(pending.hi-pending.lo)*T(0.5);
                diagnostic->stop(RootStopReason::VisitLimit,visits,0,size,
                    pending.lo,pending.hi,tolerance*max_value(T(1),abs_value(middle)));
                // 以下は停止が決定した後の診断専用評価．
                diagnostic->bounds = bound_surface(shape,o,w,pending.lo,pending.hi);
                diagnostic->has_bounds = true;
                diagnostic->f_lower = shape.implicit_value(o.at(w,pending.lo));
                diagnostic->f_upper = shape.implicit_value(o.at(w,pending.hi));
                diagnostic->has_endpoint_values = true;
            }
            return {BoundaryCode::Unresolved,T(0)};
        }
        const Segment s=stack[--size];
        const auto bounds=bound_surface(shape,o,w,s.lo,s.hi);
        if(bounds.value.excludes_zero()) continue;
        const T mid=s.lo+(s.hi-s.lo)*T(0.5);
        const T width_tolerance=tolerance*max_value(T(1),abs_value(mid));
        if(bounds.derivative.excludes_zero())
        {
            T lo=s.lo,hi=s.hi;
            T f_lo=shape.implicit_value(o.at(w,lo));
            T f_hi=shape.implicit_value(o.at(w,hi));
            if(f_lo==T(0) && lo>start) return {BoundaryCode::Hit,lo};
            if(f_hi==T(0) && hi>start) return {BoundaryCode::Hit,hi};
            if((f_lo<T(0))!=(f_hi<T(0)))
            {
                T x=mid;
                for(unsigned iteration=0;iteration<96;++iteration)
                {
                    const Vec3T<T> p=o.at(w,x);
                    const T f=shape.implicit_value(p);
                    if(f==T(0) || hi-lo<=width_tolerance) return {BoundaryCode::Hit,x};
                    if((f<T(0))==(f_lo<T(0))) {lo=x;f_lo=f;} else {hi=x;f_hi=f;}
                    const T derivative=shape.gradient(p).dot(w);
                    const T candidate=x-f/derivative;
                    const T guard=(hi-lo)*T(0.1);
                    x=candidate>lo+guard && candidate<hi-guard?candidate:lo+(hi-lo)*T(0.5);
                }
                if(diagnostic != nullptr)
                {
                    diagnostic->stop(RootStopReason::BracketIterationLimit,visits,96,size,
                        lo,hi,width_tolerance);
                    // 最後の bracket 自体の区間値を記録する．再評価は探索に使わない．
                    diagnostic->bounds = bound_surface(shape,o,w,lo,hi);
                    diagnostic->has_bounds = true;
                    diagnostic->f_lower = f_lo;
                    diagnostic->f_upper = f_hi;
                    diagnostic->has_endpoint_values = true;
                }
                return {BoundaryCode::Unresolved,T(0)};
            }
        }
        if(s.hi-s.lo<=width_tolerance || mid==s.lo || mid==s.hi || size+2>64)
        {
            if(diagnostic != nullptr)
            {
                const RootStopReason reason = s.hi-s.lo<=width_tolerance
                    ? RootStopReason::WidthLimit
                    : (mid==s.lo || mid==s.hi)
                        ? RootStopReason::MidpointStagnation : RootStopReason::StackLimit;
                diagnostic->stop(reason,visits,0,size,s.lo,s.hi,width_tolerance);
                diagnostic->bounds = bounds;
                diagnostic->has_bounds = true;
                // 追加の評価は終了が決定した後だけ．探索の分岐には利用しない．
                diagnostic->f_lower = shape.implicit_value(o.at(w,s.lo));
                diagnostic->f_upper = shape.implicit_value(o.at(w,s.hi));
                diagnostic->has_endpoint_values = true;
            }
            return {BoundaryCode::Unresolved,T(0)};
        }
        stack[size++]={mid,s.hi};
        stack[size++]={s.lo,mid};
    }
    return {BoundaryCode::Miss,T(0)};
}
}

struct RaindropIntersector
{
    const RaindropShape& shape;
    HOST_DEVICE BoundaryHit intersect(
        const BoundaryRay& ray,
        detail::BoundarySolveDiagnostics* diagnostic = nullptr) const noexcept
    {
        if(diagnostic != nullptr) *diagnostic = {};
        BoundaryHit hit{};
        const auto root=detail::solve_surface<float>(shape,ray,
            diagnostic != nullptr ? &diagnostic->fp32 : nullptr);
        if(diagnostic != nullptr) diagnostic->fp32.result = root.code;
        bool fallback=root.code==BoundaryCode::Unresolved;
        float distance=root.distance;
        if(root.code==BoundaryCode::Hit)
        {
            const Vec3 p=ray.physical_origin.at(ray.direction,distance);
            const float incidence=::fabsf(shape.gradient(p).normalized().dot(ray.direction));
            fallback=incidence<0.05f; // grazing の条件悪化だけを選択的に再評価．
        }
        hit.code=root.code;
        if(fallback)
        {
            if(diagnostic != nullptr) diagnostic->used_fp64 = true;
            const auto precise=detail::solve_surface<double>(shape,ray,
                diagnostic != nullptr ? &diagnostic->fp64 : nullptr);
            if(diagnostic != nullptr) diagnostic->fp64.result = precise.code;
            hit.code=precise.code;
            distance=static_cast<float>(precise.distance);
            hit.diagnostics|=IntersectionFallback;
            if(precise.code==BoundaryCode::Hit)
            {
                const auto p=ray.physical_origin.cast<double>().at(ray.direction.cast<double>(),precise.distance);
                hit.position=p.cast<float>();
                hit.outward_normal=shape.gradient(p).normalized().cast<float>().normalized();
            }
        }
        else if(root.code==BoundaryCode::Hit)
        {
            hit.position=ray.physical_origin.at(ray.direction,distance);
            hit.outward_normal=shape.gradient(hit.position).normalized();
        }
        hit.distance=distance;
        return hit;
    }
};
}
