#pragma once

#include <cuda.h>
#include <filesystem>

namespace rainbow
{

// MEMO: このプロジェクトでは，所有権や不変条件をもつ型は class として実装する

// CUDA Driver API と OptiX が共有する実行環境を管理する RAII クラス
// 指定された CUDA device の primary context に対する retain の参照を１つ保持し，
// その context 上に CUDA stream を一つ作成する．
// primary context 自体を排他的に所有するわけではなく，同じプロセス内の CUDA rRuntime API や OptiX と共有される可能性がある
class CudaContext final
{
public:
    explicit CudaContext(int device_ordinal = 0);       // CUDA 実行環境を構築するための device ordinal であることを呼び出し側に明示する
    ~CudaContext() noexcept;                            // デストラクタが例外を出しても外に伝播させない

    // コピーの禁止 (同じ CUDA handle を複数のオブジェクトが破棄しないようにする)
    CudaContext(const CudaContext&) = delete;           // コピーコンストラクタの禁止
    CudaContext& operator=(const CudaContext&) = delete;// コピー代入の禁止

    // move の禁止　(プログラムの上位で一度構築し，参照で渡す設計だから)
    CudaContext(CudaContext&&) = delete;                // ムーブコンストラクタの禁止
    CudaContext& operator=(CudaContext&&) = delete;     // ムーブ代入の禁止
    
    // この context を，呼び出したホストスレッドの current context にする
    void make_current() const;

    [[nodiscard]]   // 保持している CUdevice を返す observer なので，戻り値が無視されたときに警告を出してやる
    CUdevice device() const noexcept
    {
        return device_;
    }

    [[nodiscard]]
    CUcontext handle() const noexcept
    {
        return context_;
    }

    [[nodiscard]]
    CUstream stream() const noexcept
    {
        return stream_;
    }

private:
    CUdevice device_ = 0;
    CUcontext context_ = nullptr;
    CUstream stream_ = nullptr;

    bool primary_context_retained_ = false;
};

void run_cuda_driver_smoke_test(
    const CudaContext& cuda_context,
    const std::filesystem::path& fatbin_path
);
}