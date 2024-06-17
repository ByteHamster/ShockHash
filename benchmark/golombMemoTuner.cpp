#include <vector>
#include <iostream>
#include "XorShift64.h"
#include "ShockHash2-internal.h"

template <template<size_t> class T, size_t leafSize>
void dispatchLeafSize() {
    if constexpr (leafSize > 1) {
        dispatchLeafSize<T, leafSize - 1>();
    }
    if constexpr (leafSize == 1) {
        std::cout << " 1, " << std::flush; // Additional for 0
    }
    if constexpr (leafSize <= 5) {
        std::cout << " 1, " << std::flush;
        return;
    }

    std::vector<uint64_t> leaf(leafSize);
    util::XorShift64 prng;
    std::vector<size_t> seeds;
    size_t iterations = 10;
    if (leafSize < 50) {
        iterations = 100000;
    } else if (leafSize < 80) {
        iterations = 10000;
    } else if (leafSize < 130) {
        iterations = 100;
    }
    seeds.reserve(iterations);
    for (size_t i = 0; i < iterations; i++) {
        for (size_t k = 0; k < leafSize; k++) {
            leaf[k] = prng();
        }
        seeds.push_back(T<leafSize>::findSeed(leaf));
    }

    size_t spaceBest = std::numeric_limits<size_t>::max();
    size_t lowerBest = 0;
    for (size_t lower = 0; lower < 60; lower++) {
        size_t spaceBits = 0;
        for (size_t seed : seeds) {
            spaceBits += (seed >> lower) + 1 + lower;
        }
        if (spaceBits < spaceBest) {
            spaceBest = spaceBits;
            lowerBest = lower;
        }
    }
    std::cout << (lowerBest < 10 ? " " : "") << lowerBest << ", " << std::flush;
    if (leafSize % 10 == 9) {
        std::cout << "// " << (leafSize / 10) * 10 << ".." << leafSize << std::endl;
    }
}

template <size_t leafSize>
using ShockHash2 = std::conditional_t<(leafSize >= 10),
        shockhash::BijectionsShockHash2<leafSize, true, shockhash::QuadSplitCandidateFinderBuckets>,
        shockhash::BijectionsShockHash2<leafSize, true, shockhash::BasicSeedCandidateFinder>>;

int main() {
    constexpr size_t maxLeafSize = 138;
    dispatchLeafSize<ShockHash2, maxLeafSize>();
    if constexpr ((maxLeafSize % 10) != 9) {
        std::cout << "// " << (maxLeafSize / 10) * 10 << ".." << maxLeafSize << std::endl;
    }
}
