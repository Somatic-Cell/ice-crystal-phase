#pragma once

#include <filesystem>
#include <vector>

namespace rainbow
{

// fatbin / OptiX IR を対象としたバイナリファイルの読み込みを行う関数
// CUDA / OptiX の SDK には依存しない
[[nodiscard]]
std::vector<char> read_binary_file(const std::filesystem::path& path);
}