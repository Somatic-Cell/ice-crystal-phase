#include <ice_crystal/hex_trace_cuda.hpp>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

using namespace iceCrystal;
namespace
{
void check(bool b,const char* m){if(!b)throw std::runtime_error(m);}
void near(double a,double b,double t,const char* m)
{check(finite_value(a)&&finite_value(b)&&::fabs(a-b)<=t,m);}
void vector_near(Vec3 a,Vec3 b,double t,const char* m)
{check((a-b).max_abs()<=t,m);}
void compare(const TraceBatch& cpu,const TraceBatch& gpu)
{
    check(cpu.offsets==gpu.offsets,"CPU/CUDA CSR offset mismatch");
    check(cpu.audits.size()==gpu.audits.size(),"CPU/CUDA audit count");
    for(std::size_t i=0;i<cpu.audits.size();++i)
    {
        const auto& a=cpu.audits[i];const auto& b=gpu.audits[i];
        check(a.status==b.status&&a.output_count==b.output_count&&a.tir_count==b.tir_count&&a.internal_hits==b.internal_hits,"CPU/CUDA status/path topology");
        near(a.escaped_power,b.escaped_power,5e-12,"CPU/CUDA escaped power");
        near(a.unresolved_power,b.unresolved_power,5e-12,"CPU/CUDA unresolved power");
    }
    check(cpu.outgoing.size()==gpu.outgoing.size(),"CPU/CUDA outgoing count");
    for(std::size_t i=0;i<cpu.outgoing.size();++i)
    {
        const auto& a=cpu.outgoing[i];const auto& b=gpu.outgoing[i];
        check(a.incident_sample_id==b.incident_sample_id&&a.kind==b.kind&&a.entry_face==b.entry_face&&
              a.exit_face==b.exit_face&&a.internal_reflections==b.internal_reflections,"CPU/CUDA output topology");
        near(a.power_fraction,b.power_fraction,5e-12,"CPU/CUDA sample power");
        vector_near(a.direction,b.direction,5e-12,"CPU/CUDA direction");
        vector_near(a.position_mm,b.position_mm,5e-10,"CPU/CUDA position");
        vector_near(a.frame.e0,b.frame.e0,5e-12,"CPU/CUDA e0");
        vector_near(a.frame.e1,b.frame.e1,5e-12,"CPU/CUDA e1");
        near(a.internal_length_mm,b.internal_length_mm,5e-9,"CPU/CUDA path length");
        for(unsigned c=0;c<2;++c)
        {
            near(a.response.column[c].x.re,b.response.column[c].x.re,5e-12,"CPU/CUDA J real");
            near(a.response.column[c].x.im,b.response.column[c].x.im,5e-12,"CPU/CUDA J imag");
            near(a.response.column[c].y.re,b.response.column[c].y.re,5e-12,"CPU/CUDA J real");
            near(a.response.column[c].y.im,b.response.column[c].y.im,5e-12,"CPU/CUDA J imag");
        }
    }
}
}
int main(int argc,char** argv)
{
    try
    {
        const auto module=argc==2?std::filesystem::path(argv[1]):
            std::filesystem::absolute(argv[0]).parent_path()/"modules"/"hex_trace.fatbin";
        rainbow::CudaContext context(0);HexPrismTracer tracer(context);tracer.load_module(module);
        for(const auto k:std::vector<Vec3>{{0,-1,0},{1,0,0},{1,-0.4,0.3},{0.17,0.84,-0.33}})
        {
            auto s=make_trace_settings({1,2},k,1,1.31,1024,12345,256,1e-12);
            auto cpu=trace_batch_cpu(s,0,1024);tracer.trace(s,0,1024);auto gpu=tracer.download();
            require_accepted(summarize(cpu.audits));require_accepted(summarize(gpu.audits));compare(cpu,gpu);
        }
        auto limited=make_trace_settings({1,2},{0,-1,0},1,1.31,512,7,1,1e-12);
        auto cpu=trace_batch_cpu(limited,17,333);tracer.trace(limited,17,333);auto gpu=tracer.download();
        check(!summarize(gpu.audits).accepted,"CUDA truncation wrongly accepted");compare(cpu,gpu);
        auto matched=limited;matched.interior_index=1;matched.max_internal_hits=256;matched.residual_power_tolerance=0;
        cpu=trace_batch_cpu(matched,0,512);tracer.trace(matched,0,512);gpu=tracer.download();
        require_accepted(summarize(gpu.audits));compare(cpu,gpu);
        std::cout<<"PASS CUDA: six configurations, CPU comparison including Jones entries, offsets and failures\n";return 0;
    }
    catch(const std::exception& e){std::cerr<<"CUDA test FAILED (not skipped): "<<e.what()<<'\n';return 1;}
}
