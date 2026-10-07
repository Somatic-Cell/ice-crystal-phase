#include <rainbow/npy_writer.hpp>
#include <rainbow/phase_cdf_math.hpp>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void emit(const std::filesystem::path& path)
{
    using namespace rainbow;
    constexpr std::uint32_t nt=17,np=13;
    const auto edges=make_phase_u_edges(nt);
    std::vector<double> phi(np+1,0),theta(np*(nt+1),0),weights(nt*np),columns(np,0);
    for(unsigned j=0;j<np;++j)
        for(unsigned i=0;i<nt;++i)
        {
            const double d=(j==0 || j==np-1 || (i+2*j)%5==0)?0.0:0.3+(i+1.0)*(i+1.0)+0.01*j;
            weights[j*nt+i]=d*(edges[i+1]-edges[i]);columns[j]+=weights[j*nt+i];
        }
    for(unsigned j=0;j<np;++j)phi[j+1]=phi[j]+columns[j];
    const double total=phi.back();for(auto& x:phi)x/=total;phi.back()=1;
    for(unsigned j=0;j<np;++j)
    {
        const auto base=j*(nt+1u);
        for(unsigned i=0;i<nt;++i)theta[base+i+1]=theta[base+i]+weights[j*nt+i];
        for(unsigned i=0;i<=nt;++i)theta[base+i]=columns[j]>0?theta[base+i]/columns[j]:edges[i];
        theta[base]=0;theta[base+nt]=1;
    }
    DatasetDirectory directory(path);
    const std::array<std::uint64_t,1> a{np+1},e{nt+1};
    const std::array<std::uint64_t,2> b{np,nt+1};
    {NpyFloat64Writer file(directory.staging()/"phi_cdf.npy",a);file.append(phi);file.finish();}
    {NpyFloat64Writer file(directory.staging()/"theta_given_phi_cdf.npy",b);
        file.append(std::span<const double>(theta).first(7));
        file.append(std::span<const double>(theta).subspan(7));file.finish();}
    {NpyFloat64Writer file(directory.staging()/"u_edges.npy",e);file.append(edges);file.finish();}
    const std::array<std::uint64_t,2> mass_shape{np,nt};
    for(auto& x:weights)x/=total;
    {NpyFloat64Writer file(directory.staging()/"expected_mass.npy",mass_shape);file.append(weights);file.finish();}
    const std::array<std::uint64_t,8> bits{
        0,0x8000000000000000ull,1,0x3fefffffffffffffull,
        0x3ff0000000000000ull,0x7fefffffffffffffull,0x7ff0000000000000ull,0x7ff8000000001234ull};
    std::array<double,8> values{};for(unsigned i=0;i<8;++i)values[i]=std::bit_cast<double>(bits[i]);
    const std::array<std::uint64_t,1> shape{8};
    {NpyFloat64Writer file(directory.staging()/"bits.npy",shape);file.append(values);file.finish();}
    {
        std::ofstream out(directory.staging()/"metadata.json");
        out << R"JSON({"schema":"rainbow.phase_cdf.numpy.v1","complete":true,"dtype":"<f8","order":"C",
"theta_count":17,"phi_count":13,"density_measure":"solid_angle_sr",
"cell_model":"constant_density_per_spherical_cell","cdf_axis_order":["phi_cell","theta_edge"],
"coordinates":"u=(1-cos(theta))/2; v=(phi+pi)/(2*pi)",
"normalization":"integral_p_domega_equals_one","input_polarization":"unpolarized",
"sampling_frame_columns":[[1,0,0],[0,1,0],[0,0,1]],
"source":"synthetic_cpp_test_not_optical_ground_truth"})JSON";
        out.close();require(!out.fail(),"Metadata close failed");
    }
    directory.commit();
    bool rejected=false;try{DatasetDirectory again(path);}catch(const std::exception&){rejected=true;}
    require(rejected,"Writer overwrote existing dataset");
}
}
#if defined(_WIN32)
int wmain(int argc,wchar_t** argv)
#else
int main(int argc,char** argv)
#endif
{
    try
    {
        // With a path argument, retain synthetic arrays for an independent Python roundtrip.
        if(argc==2){emit(std::filesystem::path(argv[1]));return 0;}
        if(argc!=1)throw std::invalid_argument("Usage: phase_numpy_cpu [new_fixture_directory]");
        const auto nonce=std::chrono::steady_clock::now().time_since_epoch().count();
        const auto path=std::filesystem::temp_directory_path()/("rainbow_npy_test_"+std::to_string(nonce));
        emit(path);
        // All stream objects have been closed before Windows cleanup.
        std::filesystem::remove_all(path);
        std::cout<<"NPY writer and saved mass fixture passed.\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
