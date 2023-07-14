#include <chrono>
#include <iostream>
#include <XorShift64.h>
#include <tlx/cmdline_parser.hpp>
#ifdef SIMD
#include "SIMDShockHash.hpp"
template <size_t l>
using ShockHash = shockhash::SIMDShockHash<l>;
template <size_t l>
using ShockHashRotate = shockhash::SIMDShockHash<l;
#else
#include "ShockHash.h"
template <size_t l>
using ShockHash = shockhash::ShockHash<l, false>;
template <size_t l>
using ShockHashRotate = shockhash::ShockHash<l, true>;
#endif

#define DO_NOT_OPTIMIZE(value) asm volatile("" : : "r,m"(value) : "memory")

bool rotate = false;
size_t numObjects = 1e6;
size_t numQueries = 1e6;
size_t leafSize = 20;
size_t bucketSize = 2000;

template<typename HashFunc>
void construct() {
    auto time = std::chrono::system_clock::now();
    long seed = std::chrono::duration_cast<std::chrono::milliseconds>(time.time_since_epoch()).count();
    std::cout<<"Generating input data (Seed: "<<seed<<")"<<std::endl;
    util::XorShift64 prng(seed);
	std::vector<sux::function::hash128_t> keys;
    for (size_t i = 0; i < numObjects; i++) {
        keys.push_back(sux::function::hash128_t(prng(), prng()));
    }

    std::cout<<"Constructing"<<std::endl;
    auto beginConstruction = std::chrono::high_resolution_clock::now();
    HashFunc hashFunc(keys, bucketSize);
    unsigned long constructionDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - beginConstruction).count();

    std::cout<<"Testing"<<std::endl;
    std::vector<bool> taken(keys.size(), false);
    for (auto key : keys) {
        size_t hash = hashFunc(key);
        if (taken[hash]) {
            std::cerr << "Collision!" << std::endl;
            exit(1);
        } else if (hash > numObjects) {
            std::cerr << "Out of range!" << std::endl;
            exit(1);
        }
        taken[hash] = true;
    }

    std::cout<<"Querying"<<std::endl;
    uint64_t h = 0;
    auto beginQueries = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < numQueries; i++) {
        h ^= hashFunc(sux::function::hash128_t(prng(), prng() ^ h));
    }
    auto queryDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - beginQueries).count();
    DO_NOT_OPTIMIZE(h);

    hashFunc.printBits();

    std::cout << "RESULT"
              #ifdef SIMD
              << " method=SIMD"
              #else
              << " method=plain"
              #endif
                << (rotate ? "Rotate" : "")
              << " l=" << leafSize
              << " b=" << bucketSize
              << " N=" << numObjects
              << " numQueries=" << numQueries
              << " queryTimeMilliseconds=" << queryDurationMs
              << " constructionTimeMilliseconds=" << constructionDurationMs
              << " bitsPerElement=" << (double) hashFunc.getBits() / numObjects
              << std::endl;
}

template <template<size_t> class RecSplit, size_t I>
void dispatchLeafSize(size_t param) {
    if constexpr (I <= 2) {
        std::cerr<<"The parameter "<<param<<" for the leaf size was not compiled into this binary."<<std::endl;
    } else if (I == param) {
        construct<RecSplit<I>>();
    } else {
        dispatchLeafSize<RecSplit, I - 1>(param);
    }
}

int main(int argc, const char* const* argv) {
    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "numObjects", numObjects, "Number of objects to construct with");
    cmd.add_bytes('q', "numQueries", numQueries, "Number of queries to measure");
    cmd.add_bytes('l', "leafSize", leafSize, "Leaf size to construct");
    cmd.add_bytes('b', "bucketSize", bucketSize, "Bucket size to construct");
    cmd.add_bool('r', "rotate", rotate, "Apply rotation fitting");

    if (!cmd.process(argc, argv)) {
        return 1;
    }

    if (rotate) {
        dispatchLeafSize<ShockHashRotate, shockhash::MAX_LEAF_SIZE>(leafSize);
    } else {
        dispatchLeafSize<ShockHash, shockhash::MAX_LEAF_SIZE>(leafSize);
    }
    return 0;
}
