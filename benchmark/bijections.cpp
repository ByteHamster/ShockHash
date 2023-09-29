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
    T t;
    std::vector<uint64_t> leaf(leafSize);
    util::XorShift64 prng;
    double spaceSum = 0;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (size_t i = 0; i < iterations; i++) {
        for (size_t k = 0; k < leafSize; k++) {
            leaf[k] = prng();
        }
        spaceSum += t.calculateBijection(leaf);
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double usPerLeaf = 0.001 * (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / (double)iterations;
    double spaceAvg = (double)spaceSum / iterations;
    std::cout << " usPerLeaf=" << usPerLeaf
              << " avgLogHashValue=" << spaceAvg
              << " spacePerObject=" << spaceAvg / leafSize;
}

template<uint8_t leafSize>
void test(size_t iterations, size_t iterationsSlow) {
    (void) iterationsSlow;
    /*std::cout<<"RESULT l="<<(int)leafSize<<" name=RecSplit "<<std::flush;
    testSingle<leafSize, BijectionsRecSplit<leafSize>>(iterationsSlow);
    std::cout<<std::endl;

    if constexpr (leafSize <= 16) {
        std::cout<<"RESULT l="<<(int)leafSize<<" name=RecSplitRotate "<<std::flush;
        testSingle<leafSize, BijectionsRotate<leafSize, false>>(iterations);
        std::cout << std::endl;

        std::cout<<"RESULT l="<<(int)leafSize<<" name=RecSplitRotateLookup "<<std::flush;
        testSingle<leafSize, BijectionsRotate<leafSize, true>>(iterations);
        std::cout << std::endl;
    }

    if constexpr (leafSize <= 64) {
        std::cout<<"RESULT l="<<(int)leafSize<<" name=ShockHash "<<std::flush;
        testSingle<leafSize, BijectionsShockHash1<leafSize>>(iterations / 10);
        std::cout<<std::endl;

        std::cout<<"RESULT l="<<(int)leafSize<<" name=ShockHashRotate "<<std::flush;
        testSingle<leafSize, BijectionsShockHash1Rotate<leafSize>>(iterations / 10);
        std::cout<<std::endl;
    }*/

    /*std::cout<<"RESULT l="<<(int)leafSize<<" name=ShockHash2"<<std::flush;
    testSingle<leafSize, BijectionsShockHash2<leafSize, false, BasicSeedCandidateFinder>>(iterations);
    std::cout<<std::endl;*/

    /*std::cout<<"RESULT l="<<(int)leafSize<<" name=ShockHash2Filter"<<std::flush;
    testSingle<leafSize, BijectionsShockHash2<leafSize, true, BasicSeedCandidateFinder>>(iterations);
    std::cout<<std::endl;*/

    /*std::cout<<"RESULT l="<<(int)leafSize<<" name=ShockHash2RotateFilter"<<std::flush;
    testSingle<leafSize, BijectionsShockHash2<leafSize, true, RotatingSeedCandidateFinder>>(iterations);
    std::cout<<std::endl;*/

    /*std::cout<<"RESULT l="<<(int)leafSize<<" name=ShockHash2QuadSplitFilter"<<std::flush;
    testSingle<leafSize, BijectionsShockHash2<leafSize, true, QuadSplitCandidateFinderList>>(iterations);
    std::cout<<std::endl;*/

    std::cout<<"RESULT l="<<(int)leafSize<<" name=ShockHash2QuadSplitTreeFilter"<<std::flush;
    testSingle<leafSize, BijectionsShockHash2<leafSize, true, QuadSplitCandidateFinderTree>>(iterations);
    std::cout<<std::endl;

    std::cout<<std::endl;
}

int main() {
    /*test<4> (6000000, 100000);
    test<8> (6000000, 100000);
    test<12>(6000000, 6000);
    test<16>(4000000, 500);
    test<20>(1000000, 100);
    test<24>(500000, 1);
    test<28>(100000, 1);
    test<32>(100000, 1);*/
    test<36>(50000, 1);
    test<40>(10000, 0);
    test<44>(10000, 0);
    test<48>(5000, 0);
    test<52>(5000, 0);
    test<56>(5000, 0);
    test<60>(1000, 0);
    test<64>(1000, 0);
    test<70>(1000, 0);
    test<80>(1000, 0);
    test<90>(100, 0);
    test<100>(30, 0);
    test<110>(20, 0);
    test<120>(50, 0);
    test<128>(50, 0);
}