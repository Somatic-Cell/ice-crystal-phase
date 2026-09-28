#include <rainbow/read_binary_file.hpp>

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace rainbow
{

std::vector<char> read_binary_file(
    const std::filesystem::path& path
){
    // stream の寿命でファイルを閉じるため，独自のファイル所有クラスは作らない

    // ファイルの末尾を調べる
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    
    // 開けるかどうか
    if(!input)
    {
        throw std::runtime_error("Failed to open Optix IR: " + path.string());
    }

    // ファイルサイズを読めるかどうか
    const std::streampos end_position = input.tellg();
    if(end_position == std::streampos(-1))
    {
        throw std::runtime_error("Failed to read OptiX IR size: " + path.string());
    }

    // 空のファイルは拒否
    const std::streamoff file_byte_size =
        static_cast<std::streamoff>(end_position);
    if(file_byte_size <= 0)
    {
        throw std::runtime_error("OptiX IR is empty: " + path.string());
    }

    // std::uintmax_t で表現しきれないサイズも拒否
    const std::uintmax_t byte_size = static_cast<std::uintmax_t>(file_byte_size);
    if(
        byte_size > static_cast<std::uintmax_t>(std::numeric_limits<std::size_t>::max())
        || byte_size > static_cast<std::uintmax_t>(std::numeric_limits<std::streamsize>::max())
    )
    {
        throw std::runtime_error("Optix IR is too large: " + path.string());
    }

    // 先頭から byte_size 分のデータを読む
    std::vector<char> bytes(static_cast<std::size_t>(byte_size));
    input.seekg(0, std::ios::beg);
    if(!input.read(bytes.data(), static_cast<std::streamsize>(bytes.size())))
    {
        throw std::runtime_error("Failed to read OptiX IR: " + path.string());
    }

    // vector の値を返す．GPU module のロードや初期化は呼び出し側が行う．
    return bytes;

}

} // namespace rainbow