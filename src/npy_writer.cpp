#include <rainbow/npy_writer.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>

namespace rainbow
{
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
static_assert(std::endian::native == std::endian::little ||
              std::endian::native == std::endian::big);

NpyFloat64Writer::NpyFloat64Writer(const std::filesystem::path& path,
                                 std::span<const std::uint64_t> shape)
{
    if(shape.empty() || shape.size() > 8)
        throw std::invalid_argument("NPY: require 1..8 dimensions.");
    expected_ = 1;
    std::ostringstream dimensions;
    dimensions.imbue(std::locale::classic());
    dimensions << '(';
    for(const auto n : shape)
    {
        if(n == 0 || expected_ > (std::numeric_limits<std::uint64_t>::max)() / n)
            throw std::length_error("NPY: invalid/overflowing shape.");
        expected_ *= n;
        dimensions << n << ", ";
    }
    dimensions << ')';
    if(expected_ > ((std::numeric_limits<std::uint64_t>::max)() - 65546u) / 8u)
        throw std::length_error("NPY: byte size overflows.");
    std::string header = "{'descr': '<f8', 'fortran_order': False, 'shape': "
                       + dimensions.str() + ", }";
    const auto padding = (64u - ((10u + header.size() + 1u) % 64u)) % 64u;
    header.append(padding, ' ');
    header.push_back('\n');
    if(header.size() > 65535u) throw std::length_error("NPY: v1 header is too large.");
    // Directory ownership prevents concurrent writers. Do not use this writer
    // directly on an important existing path: it intentionally truncates a file.
    output_.open(path, std::ios::binary | std::ios::trunc);
    if(!output_) throw std::runtime_error("NPY: cannot open output file.");
    const std::array<unsigned char, 10> prefix{
        0x93, 'N', 'U', 'M', 'P', 'Y', 1, 0,
        static_cast<unsigned char>(header.size() & 255u),
        static_cast<unsigned char>((header.size() >> 8u) & 255u)};
    output_.write(reinterpret_cast<const char*>(prefix.data()), 10);
    output_.write(header.data(), static_cast<std::streamsize>(header.size()));
    if(!output_) throw std::runtime_error("NPY: header write failed.");
}

void NpyFloat64Writer::append(std::span<const double> values)
{
    if(finished_ || values.size() > expected_ - written_)
        throw std::logic_error("NPY: excess data or already finished.");
    constexpr std::size_t chunk_elements = 65536;
    for(std::size_t offset = 0; offset < values.size();)
    {
        const auto n = (std::min)(chunk_elements, values.size() - offset);
        if constexpr(std::endian::native == std::endian::little)
        {
            output_.write(reinterpret_cast<const char*>(values.data() + offset),
                          static_cast<std::streamsize>(n * sizeof(double)));
        }
        else
        {
            std::array<unsigned char, chunk_elements * 8> bytes{};
            for(std::size_t i = 0; i < n; ++i)
            {
                const auto bits = std::bit_cast<std::uint64_t>(values[offset + i]);
                for(unsigned b = 0; b < 8; ++b)
                    bytes[8*i+b] = static_cast<unsigned char>((bits >> (8*b)) & 255u);
            }
            output_.write(reinterpret_cast<const char*>(bytes.data()),
                          static_cast<std::streamsize>(n * 8));
        }
        if(!output_) throw std::runtime_error("NPY: data write failed (check free disk space).");
        offset += n;
        written_ += n;
    }
}
void NpyFloat64Writer::finish()
{
    if(finished_ || written_ != expected_)
        throw std::logic_error("NPY: incorrect element count or repeated finish.");
    output_.flush();
    if(!output_) throw std::runtime_error("NPY: flush failed.");
    output_.close();
    if(output_.fail()) throw std::runtime_error("NPY: close failed.");
    finished_ = true;
}

DatasetDirectory::DatasetDirectory(const std::filesystem::path& destination)
    : destination_(destination)
{
    if(destination_.empty() || destination_.filename().empty() ||
       destination_.filename() == "." || destination_.filename() == ".." ||
       destination_.extension() == ".part")
        throw std::invalid_argument("Invalid dataset destination.");
    staging_ = destination_;
    staging_ += ".part";
    if(std::filesystem::exists(destination_))
        throw std::runtime_error("Dataset already exists; refusing to overwrite.");
    if(!destination_.parent_path().empty())
        std::filesystem::create_directories(destination_.parent_path());
    if(!std::filesystem::create_directory(staging_))
        throw std::runtime_error("Staging directory exists; inspect/retire the previous .part directory.");
}
void DatasetDirectory::commit()
{
    if(committed_ || !std::filesystem::is_regular_file(staging_ / "metadata.json"))
        throw std::logic_error("Dataset has no complete metadata or was already committed.");
    if(std::filesystem::exists(destination_))
        throw std::runtime_error("Destination appeared while writing; refusing to replace it.");
    std::filesystem::rename(staging_, destination_);
    committed_ = true;
}
} // namespace rainbow
