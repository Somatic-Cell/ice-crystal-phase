#pragma once

#include <rainbow/bilinear_patch_intersection.hpp>
#include <cstdint>
#include <cstddef>
#include <type_traits>

namespace rainbow
{

// 幾何の問い合わせ結果．強度・位相関数ではない．
enum PatchQueryFlags : std::uint32_t
{
    PatchQueryNone = 0,
    PatchQueryUnresolved = 1u << 0,       // 少なくとも一候補の孤立解を確定できない．
    PatchQueryInvalidData = 1u << 1,
    PatchQueryCountMismatch = 1u << 2,    // count/fill 不一致．容量を越えた書込みは行わない．
    PatchQueryRefinementHit = 1u << 3,    // NeedsRefinement に少なくとも一つ交差．
    PatchQueryBoundaryHit = 1u << 4,      // 共有辺/corner の光学 ownership は未解決．
    PatchQueryFp64Used = 1u << 5,
    PatchQueryDuplicateReport = 1u << 6,  // 同じ primitive/root の重複を除去した．
    PatchQuerySingularHit = 1u << 7,
    PatchQueryCountOverflow = 1u << 8
};
inline constexpr std::uint32_t patch_query_error_mask =
    PatchQueryUnresolved | PatchQueryInvalidData | PatchQueryCountMismatch | PatchQueryCountOverflow;

struct PatchQueryHit
{
    std::uint32_t compact_index = 0;
    std::uint32_t patch_id = 0;     // 元の family/cell 番号（compact index と区別）．
    std::uint32_t root_index = 0;   // 同じ patch 内の u,v 辞書順．二つの根を保持．
    std::uint32_t flags = 0;        // BilinearHitFlags．
    float u = 0.0f, v = 0.0f, t = 0.0f, residual = 0.0f;
};

struct PatchQuerySummary
{
    std::uint32_t candidate_count = 0;   // count pass が確保した容量（方向ごと）．
    std::uint32_t hit_count = 0;         // 同じ primitive/root の重複だけを除いた個数．
    std::uint32_t regular_hits = 0;
    std::uint32_t refinement_hits = 0;
    std::uint32_t boundary_hits = 0;
    std::uint32_t flags = 0;
    std::uint32_t first_problem_patch_id = 0xffffffffu;
    std::uint32_t reserved = 0;
};

// 1 direction を一つの raygen invocation だけが所有する．atomic/固定長 hit list は不要．
// count/fill 両方から同じ consume() を使い，分岐条件を別実装にしない．
struct PatchHitCollector
{
    PatchQueryHit* storage = nullptr; // count pass では nullptr．
    std::uint32_t capacity = 0;
    std::uint32_t count = 0;
    std::uint32_t flags = 0;
    bool writing = false;
    std::uint32_t first_problem_patch_id = 0xffffffffu;

    HOST_DEVICE void note(const std::uint32_t geometry_flags, const std::uint32_t patch_id = 0xffffffffu) noexcept
    {
        if((geometry_flags & (BilinearUnresolved | BilinearInvalidInput)) != 0
           && patch_id < first_problem_patch_id) first_problem_patch_id=patch_id;
        if(geometry_flags & BilinearUnresolved) flags |= PatchQueryUnresolved;
        if(geometry_flags & BilinearInvalidInput) flags |= PatchQueryInvalidData;
        if(geometry_flags & BilinearFp64Evaluation) flags |= PatchQueryFp64Used;
    }

    HOST_DEVICE void consume(
        const OutgoingPatch& patch, const std::uint32_t compact_index,
        const BilinearRayHits& intersections) noexcept
    {
        note(intersections.flags,patch.patch_id);
        for(std::uint32_t root=0;root<intersections.count;++root)
        {
            const auto& h=intersections.hits[root];
            if(count==0xffffffffu) {flags|=PatchQueryCountOverflow;return;}
            if(writing)
            {
                if(count>=capacity || storage==nullptr) flags|=PatchQueryCountMismatch;
                else storage[count]={compact_index,patch.patch_id,root,h.flags,h.u,h.v,h.t,h.residual};
            }
            ++count;
        }
    }

    // ソートは「近い順」ではなく patch_id/root_index 順．BVH の訪問順序に依存しない．
    // 通常の小さな交差集合向けの in-place insertion sort．固定個数で打切らない．
    [[nodiscard]] HOST_DEVICE PatchQuerySummary finish(
        const OutgoingPatch* patches, const std::uint32_t patch_count) noexcept
    {
        PatchQuerySummary result{};
        result.candidate_count=capacity;
        if(count!=capacity) flags|=PatchQueryCountMismatch;
        result.flags=flags;
        result.first_problem_patch_id=first_problem_patch_id;
        if(count>capacity || (capacity!=0 && storage==nullptr)) return result;
        for(std::uint32_t i=1;i<count;++i)
        {
            const auto value=storage[i]; std::uint32_t j=i;
            while(j!=0 && less(value,storage[j-1])) {storage[j]=storage[j-1];--j;}
            storage[j]=value;
        }
        std::uint32_t kept=0;
        for(std::uint32_t i=0;i<count;++i)
        {
            const auto h=storage[i];
            if(kept!=0 && h.compact_index==storage[kept-1].compact_index
                       && h.root_index==storage[kept-1].root_index)
            {result.flags|=PatchQueryDuplicateReport;continue;}
            if(h.compact_index>=patch_count) {result.flags|=PatchQueryInvalidData;continue;}
            storage[kept++]=h;
            const auto status=patches[h.compact_index].status;
            if(status==PatchCellStatus::NeedsRefinement)
            {++result.refinement_hits;result.flags|=PatchQueryRefinementHit;}
            else if(status==PatchCellStatus::RegularPositive || status==PatchCellStatus::RegularNegative)
                ++result.regular_hits;
            else result.flags|=PatchQueryInvalidData;
            if(h.flags&BilinearBoundary) {++result.boundary_hits;result.flags|=PatchQueryBoundaryHit;}
            if(h.flags&BilinearSingular) result.flags|=PatchQuerySingularHit;
        }
        result.hit_count=kept;
        return result;
    }

private:
    HOST_DEVICE static bool less(const PatchQueryHit& a,const PatchQueryHit& b) noexcept
    {
        if(a.patch_id!=b.patch_id) return a.patch_id<b.patch_id;
        if(a.compact_index!=b.compact_index) return a.compact_index<b.compact_index;
        return a.root_index<b.root_index;
    }
};

static_assert(sizeof(PatchQueryHit)==32 && sizeof(PatchQuerySummary)==32);
static_assert(std::is_trivially_copyable_v<PatchQueryHit>);
static_assert(std::is_trivially_copyable_v<PatchQuerySummary>);
} // namespace rainbow
