#pragma once
#include <ice_crystal/phase_math.hpp>
#include <array>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace iceCrystal
{
struct OrientationNode
{
    Rotation rotation;
    std::array<double,4> quaternion{};
    double weight=0; // number probability mass, NOT PDF per radian or encounter probability
    std::uint64_t samples=0,seed=0;
};
struct OrientationPlan
{
    std::vector<OrientationNode> nodes;
    std::uint64_t total_samples=0;
    double input_weight_sum=0;
    std::string source_sha256;
};
std::string sha256(std::string_view bytes);
OrientationPlan parse_orientation_plan(std::string_view csv);
OrientationPlan read_orientation_plan(const std::filesystem::path&);
} // namespace iceCrystal
