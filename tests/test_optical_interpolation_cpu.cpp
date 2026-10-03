#include <rainbow/phase_interpolation.hpp>
#include <rainbow/polarization_transport.hpp>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace
{

using rainbow::Complex32;
using rainbow::Field32;
using rainbow::PhaseCycles;
using rainbow::PhaseInterpolation;
using rainbow::PolarizationTransport;
using rainbow::TransverseFrame;
using rainbow::Vec3;

[[nodiscard]] bool close(const float a,const float b,const float tolerance=2.5e-6f)
{
    return std::fabs(a-b)<=tolerance*(1.0f+std::fmax(std::fabs(a),std::fabs(b)));
}
void require(const bool condition,const char* message)
{
    if(!condition) throw std::runtime_error(message);
}
void require_complex(const Complex32 a,const Complex32 b,const char* message)
{
    require(close(a.real,b.real)&&close(a.imag,b.imag),message);
}
void require_field(const Field32 a,const Field32 b,const char* message)
{
    require_complex(a.x,b.x,message);
    require_complex(a.y,b.y,message);
}

[[nodiscard]] long double cycles(const PhaseCycles p)
{
    return static_cast<long double>(p.turns)+static_cast<long double>(p.fraction);
}

void test_phase_corners()
{
    const PhaseCycles q00{123456,0.125f},q10{123460,0.25f};
    const PhaseCycles q01{123452,0.75f},q11{123470,0.875f};
    PhaseCycles out{};
    require(PhaseInterpolation::try_bilinear(q00,q10,q01,q11,0,0,out)
            && out.turns==q00.turns && out.fraction==q00.fraction,"phase corner 00");
    require(PhaseInterpolation::try_bilinear(q00,q10,q01,q11,1,0,out)
            && out.turns==q10.turns && out.fraction==q10.fraction,"phase corner 10");
    require(PhaseInterpolation::try_bilinear(q00,q10,q01,q11,0,1,out)
            && out.turns==q01.turns && out.fraction==q01.fraction,"phase corner 01");
    require(PhaseInterpolation::try_bilinear(q00,q10,q01,q11,1,1,out)
            && out.turns==q11.turns && out.fraction==q11.fraction,"phase corner 11");
}

void test_phase_against_long_double()
{
    const PhaseCycles q00{19127,0.12345678f};
    const PhaseCycles q10{19131,0.87654322f};
    const PhaseCycles q01{19122,0.33333334f};
    const PhaseCycles q11{19139,0.77777779f};
    const float u=0.37109375f,v=0.68359375f;
    PhaseCycles out{};
    require(PhaseInterpolation::try_bilinear(q00,q10,q01,q11,u,v,out),"phase bilinear failed");

    const long double U=u,V=v;
    const long double reference=(1-U)*(1-V)*cycles(q00)+U*(1-V)*cycles(q10)
        +(1-U)*V*cycles(q01)+U*V*cycles(q11);
    const long double error=std::fabs(cycles(out)-reference);
    require(error<2.0e-6L,"phase bilinear differs from long-double reference");

    Complex32 phasor{};
    require(out.try_unit_phasor(phasor),"phase phasor failed");
    const long double fraction=reference-std::floor(reference);
    constexpr long double two_pi=6.283185307179586476925286766559005768L;
    require(close(phasor.real,static_cast<float>(std::cos(two_pi*fraction)),4.0e-6f)
            && close(phasor.imag,static_cast<float>(std::sin(two_pi*fraction)),4.0e-6f),
            "phase phasor differs from long-double reference");
}

void test_phase_integer_shift_invariance()
{
    const PhaseCycles q00{12000,0.11f},q10{12003,0.81f};
    const PhaseCycles q01{11998,0.29f},q11{12007,0.67f};
    constexpr std::int32_t shift=1000000;
    const PhaseCycles s00{q00.turns+shift,q00.fraction};
    const PhaseCycles s10{q10.turns+shift,q10.fraction};
    const PhaseCycles s01{q01.turns+shift,q01.fraction};
    const PhaseCycles s11{q11.turns+shift,q11.fraction};
    PhaseCycles a{},b{};
    require(PhaseInterpolation::try_bilinear(q00,q10,q01,q11,0.37f,0.61f,a),"base phase failed");
    require(PhaseInterpolation::try_bilinear(s00,s10,s01,s11,0.37f,0.61f,b),"shifted phase failed");
    require(b.turns-a.turns==shift,"integer path shift not preserved");
    require(a.fraction==b.fraction,"integer path shift changed fractional phase");
}

void test_phase_invalid_parameters()
{
    const PhaseCycles q{3,0.5f};
    PhaseCycles out{9,0.25f};
    require(!PhaseInterpolation::try_bilinear(q,q,q,q,-0.01f,0.5f,out),"negative u accepted");
    require(out.turns==9&&out.fraction==0.25f,"failed interpolation modified output");
}

void test_polarization_identity()
{
    const Vec3 w{0,0,1};
    const TransverseFrame frame{{1,0,0},{0,1,0}};
    const Field32 field{{0.7f,-0.2f},{-0.1f,0.4f}};
    Field32 out{};
    require(PolarizationTransport::try_transport(w,frame,field,w,frame,out),"identity transport failed");
    require_field(out,field,"identity transport changed field");
}

void test_polarization_known_rotation()
{
    // Minimum rotation z -> x is +90 degrees around +y.  It maps +x -> -z
    // and leaves +y unchanged, so the Jones components remain unchanged in
    // the explicitly rotated target frame (-z,+y).
    const Vec3 source_w{0,0,1};
    const Vec3 target_w{1,0,0};
    const TransverseFrame source_frame{{1,0,0},{0,1,0}};
    const TransverseFrame target_frame{{0,0,-1},{0,1,0}};
    const Field32 field{{0.3f,0.8f},{-0.6f,0.2f}};
    Field32 out{};
    require(PolarizationTransport::try_transport(
        source_w,source_frame,field,target_w,target_frame,out),"known rotation failed");
    require_field(out,field,"known minimum rotation is incorrect");
}

void test_polarization_gauge_invariance()
{
    const Vec3 source_w=Vec3{0.2f,-0.3f,0.9327379f}.normalized();
    const Vec3 target_w=Vec3{-0.15f,0.25f,0.9565563f}.normalized();
    const TransverseFrame source_frame=TransverseFrame::from_direction(source_w);
    const TransverseFrame target_frame=TransverseFrame::from_direction(target_w);
    const Field32 field{{0.71f,-0.13f},{-0.22f,0.51f}};

    constexpr float angle=0.731f;
    const float c=std::cos(angle),s=std::sin(angle);
    const TransverseFrame rotated_frame{
        source_frame.e0*c+source_frame.e1*s,
        source_frame.e0*(-s)+source_frame.e1*c
    };
    const Field32 rotated_components=field.rotated_basis(c,s);

    Field32 a{},b{};
    require(PolarizationTransport::try_transport(
        source_w,source_frame,field,target_w,target_frame,a),"base gauge transport failed");
    require(PolarizationTransport::try_transport(
        source_w,rotated_frame,rotated_components,target_w,target_frame,b),"rotated gauge transport failed");
    require_field(a,b,"source-frame gauge changed transported physical field");
    require(close(a.squared_norm(),field.squared_norm(),5.0e-6f),"transport did not preserve field norm");
}

void test_scattering_frame()
{
    const Vec3 incident{1,0,0};
    const Vec3 incident_reference{0,0,1};
    const Vec3 outgoing=Vec3{-0.7f,0.2f,0.6855655f}.normalized();
    TransverseFrame frame{};
    require(PolarizationTransport::try_make_scattering_frame(
        incident,incident_reference,outgoing,frame),"scattering frame construction failed");
    require(close(frame.e0.dot(outgoing),0.0f,2.0e-6f),"scattering e0 is not transverse");
    require(close(frame.e1.dot(outgoing),0.0f,2.0e-6f),"scattering e1 is not transverse");
    require(close(frame.e0.dot(frame.e1),0.0f,2.0e-6f),"scattering frame is not orthogonal");
    require(close(frame.e0.cross(frame.e1).dot(outgoing),1.0f,3.0e-6f),"scattering frame handedness is wrong");

    const Vec3 expected=incident.cross(outgoing).normalized();
    require(close(frame.e0.dot(expected),1.0f,3.0e-6f),"scattering perpendicular axis is wrong");
}

}

int main()
{
    try
    {
        test_phase_corners();
        test_phase_against_long_double();
        test_phase_integer_shift_invariance();
        test_phase_invalid_parameters();
        test_polarization_identity();
        test_polarization_known_rotation();
        test_polarization_gauge_invariance();
        test_scattering_frame();
        std::cout<<"Optical interpolation CPU tests: passed\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e)
    {
        std::cerr<<"Error: "<<e.what()<<'\n';
        return EXIT_FAILURE;
    }
}
