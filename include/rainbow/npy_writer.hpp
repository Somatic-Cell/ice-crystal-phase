#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>

namespace rainbow
{
// Narrow NPY 1.0 writer: C-contiguous, IEEE binary64, little endian, no pickle.
// append() streams pieces; finish() checks the exact element count and close.
class NpyFloat64Writer final
{
public:
    NpyFloat64Writer(const std::filesystem::path& path,
                     std::span<const std::uint64_t> shape);
    ~NpyFloat64Writer() = default;
    NpyFloat64Writer(const NpyFloat64Writer&) = delete;
    NpyFloat64Writer& operator=(const NpyFloat64Writer&) = delete;
    void append(std::span<const double> values);
    void finish();
private:
    std::ofstream output_;
    std::uint64_t expected_ = 0, written_ = 0;
    bool finished_ = false;
};

// Creates <destination>.part exclusively. Only commit() publishes the directory.
// Failed/incomplete directories are intentionally retained, never read as records.
// One writer per destination; does not promise power-loss durability.
class DatasetDirectory final
{
public:
    explicit DatasetDirectory(const std::filesystem::path& destination);
    DatasetDirectory(const DatasetDirectory&) = delete;
    DatasetDirectory& operator=(const DatasetDirectory&) = delete;
    [[nodiscard]] const std::filesystem::path& staging() const noexcept { return staging_; }
    void commit();
private:
    std::filesystem::path destination_, staging_;
    bool committed_ = false;
};
} // namespace rainbow
