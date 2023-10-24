#include <vector>
#include <chrono>
#include <iostream>
#include "XorShift64.h"

#include "bijections/RecSplit.h"
#include "bijections/RecSplitRotate.h"
#include "bijections/ShockHash1.h"
#include "bijections/ShockHash1Rotate.h"
#include "ShockHash2-internal.h"

template<size_t leafSize, class T>
void testSingle(size_t iterations) {
    std::cout << "RESULT"
              << " l=" << (int)leafSize
              << " name=" << T::name()
              << " iterations=" << iterations
              << std::flush;
    T t;
    std::vector<uint64_t> leaf(leafSize);
    util::XorShift64 prng;
    double spaceSum = 0;
    double seedSum = 0;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (size_t i = 0; i < iterations; i++) {
        for (size_t k = 0; k < leafSize; k++) {
            leaf[k] = prng();
        }
        size_t seed = t.calculateBijection(leaf);
        seedSum += seed;
        spaceSum += ceil(log2(seed + 1));
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double usPerLeaf = 0.001 * (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / (double)iterations;
    double spaceAvg = spaceSum / iterations;
    std::cout << " usPerLeaf=" << usPerLeaf
              << " avgLogHashValue=" << spaceAvg
              << " avgHashValue=" << (seedSum / iterations)
              << " spacePerObject=" << spaceAvg / leafSize
              << std::endl;
}

template<uint8_t leafSize>
void test() {
    using namespace shockhash;
    size_t iterationsShockHash = 4 * std::pow(1.11, 128.0 - leafSize);
    //size_t iterationsRecSplit  = 4 * std::pow(1.55, 28.0 - leafSize);

    if constexpr (leafSize <= 24) {
        //testSingle<leafSize, BijectionsRecSplit<leafSize>>(iterationsRecSplit / 10);
        //testSingle<leafSize, BijectionsRotate<leafSize, false>>(iterationsRecSplit);
        //testSingle<leafSize, BijectionsRotate<leafSize, true>>(iterationsRecSplit);
    }

    if constexpr (leafSize < 64) {
        //testSingle<leafSize, BijectionsShockHash1<leafSize>>(iterationsShockHash / 10);
        //testSingle<leafSize, BijectionsShockHash1Rotate<leafSize>>(iterationsShockHash / 10);
    }

    //testSingle<leafSize, BijectionsShockHash2<leafSize, false, BasicSeedCandidateFinder>>(iterationsShockHash);
    //testSingle<leafSize, BijectionsShockHash2<leafSize, true, BasicSeedCandidateFinder>>(iterationsShockHash);
    //testSingle<leafSize, BijectionsShockHash2<leafSize, true, RotatingSeedCandidateFinder>>(iterationsShockHash);
    testSingle<leafSize, BijectionsShockHash2<leafSize, true, QuadSplitCandidateFinderList>>(iterationsShockHash);
    //testSingle<leafSize, BijectionsShockHash2<leafSize, true, QuadSplitCandidateFinderBuckets>>(iterationsShockHash);

    std::cout<<std::endl;
}

template <size_t I>
void dispatchLeafSize() {
    if constexpr (I >= 12) {
        dispatchLeafSize<I - 2>();
        test<I>();
    }
}

int main() {
    dispatchLeafSize<128>();
}