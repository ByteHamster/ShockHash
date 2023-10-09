#include "ShockHash2-precompiled.h"
#include "ShockHash2-internal.h"

namespace shockhash {

template<size_t I>
inline size_t dispatchLeafSizeQ(size_t param, size_t seed, uint64_t key, size_t retrieved) {
    if constexpr (I <= 1) {
        return 0;
    } else if (I == param) {
        using SH = std::conditional_t<(I >= 10),
                BijectionsShockHash2<I, true, QuadSplitCandidateFinderList>,
                BijectionsShockHash2<I, true, BasicSeedCandidateFinder>>;
        return SH::hash(seed, key, retrieved);
    } else {
        return dispatchLeafSizeQ<I - 1>(param, seed, key, retrieved);
    }
}

size_t shockhash2query(size_t param, size_t seed, uint64_t key, size_t retrieved) {
   return dispatchLeafSizeQ<128>(param, seed, key, retrieved);
}

template<size_t I>
inline size_t dispatchLeafSize(size_t param, std::vector<uint64_t> &leafKeys,
                               std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput) {
    if constexpr (I <= 1) {
        return 0;
    } else if (I == param) {
        using SH = std::conditional_t<(I >= 10),
                BijectionsShockHash2<I, true, QuadSplitCandidateFinderList>,
                BijectionsShockHash2<I, true, BasicSeedCandidateFinder>>;
        size_t x = SH::findSeed(leafKeys);
        SH::constructRetrieval(leafKeys, x, ribbonInput);
#ifndef NDEBUG
        SH::verify(x, leafKeys, ribbonInput);
#endif
        return x;
    } else {
        return dispatchLeafSize<I - 1>(param, leafKeys, ribbonInput);
    }
}

size_t shockhash2construct(size_t param, std::vector<uint64_t> &leafKeys,
                        std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput) {
    return dispatchLeafSize<128>(param, leafKeys, ribbonInput);
}

}
