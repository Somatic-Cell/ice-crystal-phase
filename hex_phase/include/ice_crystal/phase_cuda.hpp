#pragma once
#include <ice_crystal/phase_table.hpp>
#include <ice_crystal/hex_trace_cuda.hpp>

namespace iceCrystal
{
// GPU scatter-add is binary64, but its order is not deterministic. Validation
// compares integrated mass against independently accumulated ray audits.
class PhaseMassAccumulatorCuda final
{
public:
    PhaseMassAccumulatorCuda(const rainbow::CudaContext&,PhaseGrid);
    ~PhaseMassAccumulatorCuda() noexcept;
    PhaseMassAccumulatorCuda(const PhaseMassAccumulatorCuda&)=delete;
    PhaseMassAccumulatorCuda& operator=(const PhaseMassAccumulatorCuda&)=delete;
    void load_module(const std::filesystem::path&);
    void add(const rainbow::DeviceBuffer<OutgoingSample>&,const Rotation&,double area_weight_per_ray);
    void finish();
    [[nodiscard]] HistogramStats statistics() const;
    [[nodiscard]] HistogramMoments moments() const;
    [[nodiscard]] const rainbow::DeviceBuffer<double>& density() const;
    [[nodiscard]] std::vector<double> download_density() const; // tests / diagnostics only
private:
    void sync();
    void launch(const char*,void**,unsigned,unsigned=256);
    PhaseGrid grid_;
    const rainbow::CudaContext& context_;
    rainbow::CudaModule module_;
    rainbow::DeviceBuffer<double> edges_,mass_,density_;
    rainbow::DeviceBuffer<HistogramStats> statistics_;
    rainbow::DeviceBuffer<HistogramMoments> partials_;
    HistogramStats host_statistics_{};
    HistogramMoments host_moments_{};
    bool pending_=false,ready_=false,finished_=false;
};
} // namespace iceCrystal
