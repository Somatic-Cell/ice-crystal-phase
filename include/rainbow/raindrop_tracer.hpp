#pragma once

#include <rainbow/optix_context.hpp>
#include <rainbow/device_buffer.hpp>
#include <rainbow/raindrop_settings.hpp>
#include <cstddef>
#include <filesystem>
#include <vector>
#include <span>


namespace rainbow
{

#if defined(_MSC_VER)
#pragma warning(push)
// SBT record を直接保持するため，クラス末尾のパディングは意図したもの．
#pragma warning(disable : 4324)
#endif

    // 非球形雨粒の path-major 出射 vertex を GPU 上に所有する実処理．テスト専用型ではない．
// create_pipeline() は一度だけ，trace() は条件ごとに反復できる．
// trace() は完了待ちして戻る初期 API．後段は vertices().address() を使えば CPU 転送不要．
class RaindropTracer final
{
public:
    explicit RaindropTracer(const CudaContext& cuda_context);
    ~RaindropTracer() noexcept;
    RaindropTracer(const RaindropTracer&)=delete;
    RaindropTracer& operator=(const RaindropTracer&)=delete;
    RaindropTracer(RaindropTracer&&)=delete;
    RaindropTracer& operator=(RaindropTracer&&)=delete;

    void create_pipeline(const std::filesystem::path& optixir_path);
    // Explicit legacy single-input API. The application uses trace_unpolarized.
    void trace(const RaindropSettings& settings);
    // One geometric trace, two unit-input Jones response columns.
    // settings.incident_field belongs only to trace(); it is not a source here.
    void trace_unpolarized(const RaindropSettings& settings);
    [[nodiscard]] std::vector<Field32> download_second_input_fields() const;
    [[nodiscard]] bool is_unpolarized() const noexcept { return unpolarized_; }
    [[nodiscard]] bool has_result() const noexcept { return has_result_; }
    [[nodiscard]] const DeviceBuffer<Field32>& second_input_fields() const noexcept { return second_input_fields_; }
    [[nodiscard]] std::vector<OutgoingVertex> download_vertices() const;
    void write_csv(const std::filesystem::path& output_path) const;
    // 既に readback したデータを使う場合，巨大な vertex buffer を二重にコピーしない．
    void write_csv(const std::filesystem::path& output_path, std::span<const OutgoingVertex> output) const;
    [[nodiscard]] const DeviceBuffer<OutgoingVertex>& vertices() const noexcept {return vertices_;}
    [[nodiscard]] const RaindropTraceConfig& config() const noexcept {return config_;}
    [[nodiscard]] const OptixContext& optix_context() const noexcept
    {
        return optix_context_;
    }
    [[nodiscard]] bool close_noexcept() noexcept;
    void close();
private:
    struct alignas(OPTIX_SBT_RECORD_ALIGNMENT) HeaderRecord {char header[OPTIX_SBT_RECORD_HEADER_SIZE];};
    static_assert(sizeof(HeaderRecord)%OPTIX_SBT_RECORD_ALIGNMENT==0);
    void trace_impl(const RaindropSettings&, bool unpolarized);
    void build_drop_gas();
    void synchronize();

    const CudaContext& cuda_context_;
    OptixContext optix_context_;
    OptixModule module_=nullptr;
    OptixProgramGroup raygen_=nullptr,miss_=nullptr,hitgroup_=nullptr;
    OptixPipeline pipeline_=nullptr;
    DeviceBuffer<HeaderRecord> raygen_record_,miss_record_,hit_record_;
    DeviceBuffer<OptixAabb> aabb_;
    DeviceBuffer<std::byte> gas_scratch_,gas_output_;
    DeviceBuffer<RaindropTraceParams> device_params_;
    DeviceBuffer<RaindropTraceConfig> device_config_;
    DeviceBuffer<OutgoingVertex> vertices_;
    DeviceBuffer<Field32> second_input_fields_;
    bool unpolarized_ = false;
    // async 転送の入力は，例外経路で同期が終わるまで生存するメンバに置く．
    HeaderRecord host_raygen_{},host_miss_{},host_hit_{};
    OptixAabb host_aabb_{};
    RaindropTraceParams params_{};
    RaindropTraceConfig config_{};
    bool has_pending_work_=false,has_result_=false,is_closed_=false;
};

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

}
