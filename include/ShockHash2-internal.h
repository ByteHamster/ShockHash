#pragma once

#include <cmath>
#include <array>
#include <vector>
#include <variant>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <TinyBinaryCuckooHashTable.h>
#include <UnionFind.h>
#include "PairingFunction.h"
#include "CuckooUnionFind.h"

template <size_t leafSize, bool filter = false>
struct SeedCache {
    [[no_unique_address]]
    std::conditional_t<filter, __uint128_t, std::monostate> isolatedVertices;
    uint64_t seed;
    uint8_t hashes[leafSize];
};

template <size_t leafSize>
inline void calculateIsolatedVertices(SeedCache<leafSize, true> &seedCache) {
    seedCache.isolatedVertices = 0;
    uint64_t hitCount[leafSize] = {0};
    for (size_t i = 0; i < leafSize; i++) {
        hitCount[seedCache.hashes[i]]++;
    }
    for (size_t i = 0; i < leafSize; i++) {
        if (hitCount[seedCache.hashes[i]] == 1) {
            seedCache.isolatedVertices |= __uint128_t(1) << i;
        }
    }

    /*__uint128_t hitCount1 = 0;
    __uint128_t hitCountMoreThan1 = 0;
    for (size_t i = 0; i < leafSize; i++) {
        // hitCount[seedCache.hashes[i]]++;
        __uint128_t pos = __uint128_t(1) << seedCache.hashes[i];
        hitCountMoreThan1 |= (hitCount1 & pos);
        hitCount1 |= pos;
    }
    __uint128_t isolated = 0;
    for (size_t i = 0; i < leafSize; i++) {
        //if (hitCount[seedCache.hashes[i]] == 1) {
        //    seedCache.isolatedVertices |= __uint128_t(1) << i;
        //}
        __uint128_t pos = __uint128_t(1) << seedCache.hashes[i];
        if ((((~hitCount1) | hitCountMoreThan1) & pos) == 0) {
            isolated |= __uint128_t(1) << i;
        }
    }
    seedCache.isolatedVertices = isolated;*/
}

template <size_t leafSize, bool filter = false>
class BasicSeedCandidateFinder {
    private:
        const std::vector<uint64_t> &keys;
        size_t currentSeed = 0;
    public:
        explicit BasicSeedCandidateFinder(const std::vector<uint64_t> &keys) : keys(keys) {

        }

        inline SeedCache<leafSize, filter> next() {
            constexpr uint64_t mask = (leafSize == 128) ? ~0ul : (1ul << ((leafSize + 1) / 2)) - 1;
            SeedCache<leafSize, filter> seedCache; // NOLINT(cppcoreguidelines-pro-type-member-init)
            while (true) {
                uint64_t taken = 0;
                for (size_t i = 0; i < leafSize; i++) {
                    uint64_t hash = util::remix(keys.at(i) + currentSeed);
                    seedCache.hashes[i] = util::fastrange64(hash, (leafSize + 1) / 2);
                    taken |= 1ul << seedCache.hashes[i];
                }
                if (taken == mask) {
                    // Found a new seed candidate
                    seedCache.seed = currentSeed;
                    if constexpr (filter) {
                        calculateIsolatedVertices(seedCache);
                    }
                    currentSeed++;
                    return seedCache;
                }
                currentSeed++;
            }
        }

        static size_t hash(uint64_t key, uint64_t seed) {
            return util::fastrange64(util::remix(key + seed), (leafSize + 1) / 2);
        }
};

template <size_t leafSize, bool filter = false>
class RotatingSeedCandidateFinder {
    private:
        std::array<uint64_t, leafSize> keys = { 0 };
        uint64_t sizeSetA = 0;
        size_t currentSeed = -1;
        size_t currentRotation = (leafSize + 1) / 2;
        uint64_t takenA = 0;
        uint64_t takenB = 0;
        SeedCache<leafSize, filter> seedCache; // NOLINT(cppcoreguidelines-pro-type-member-init)
    public:
        explicit RotatingSeedCandidateFinder(const std::vector<uint64_t> &keysIn) {
            for (size_t i = 0; i < leafSize; i++) {
                if ((keysIn[i] & 1) == 0) {
                    keys.at(sizeSetA) = keysIn[i];
                    sizeSetA++;
                } else {
                    keys.at(leafSize - i + sizeSetA - 1) = keysIn[i];
                }
            }
        }

        static constexpr uint64_t rotate(uint64_t val, uint32_t x) {
            constexpr uint64_t l = (leafSize + 1) / 2;
            return ((val << x) | (val >> (l - x))) & ((1ul << l) - 1);
        }

        inline SeedCache<leafSize, filter> next() {
            constexpr uint64_t mask = (leafSize >= 127) ? ~0ul : (1ul << ((leafSize + 1) / 2)) - 1;
            while (true) {
                while (currentRotation < (leafSize + 1) / 2) {
                    if ((takenA | rotate(takenB, currentRotation)) == mask) {
                        // Found a new seed candidate
                        SeedCache<leafSize, filter> rotated = seedCache;
                        rotated.seed = currentSeed * ((leafSize + 1) / 2) + currentRotation;
                        for (size_t i = sizeSetA; i < leafSize; i++) {
                            rotated.hashes[i] = (rotated.hashes[i] + currentRotation) % ((leafSize + 1) / 2);
                        }
                        if constexpr (filter) {
                            calculateIsolatedVertices(rotated);
                        }
                        currentRotation++;
                        return rotated;
                    }
                    currentRotation++;
                }
                currentRotation = 0;
                currentSeed++;

                takenA = 0;
                takenB = 0;
                for (size_t i = 0; i < sizeSetA; i++) {
                    uint64_t hash = util::remix(keys[i] + currentSeed);
                    seedCache.hashes[i] = util::fastrange64(hash, (leafSize + 1) / 2);
                    takenA |= 1ul << seedCache.hashes[i];
                }
                for (size_t i = sizeSetA; i < leafSize; i++) {
                    uint64_t hash = util::remix(keys[i] + currentSeed);
                    seedCache.hashes[i] = util::fastrange64(hash, (leafSize + 1) / 2);
                    takenB |= 1ul << seedCache.hashes[i];
                }
            }
        }

        static size_t hash(uint64_t key, uint64_t seed) {
            size_t hashSeed = seed / ((leafSize + 1) / 2);
            size_t rotation = seed % ((leafSize + 1) / 2);
            size_t baseHash = util::fastrange64(util::remix(key + hashSeed), (leafSize + 1) / 2);
            if ((key & 1) == 0) {
                return baseHash;
            } else {
                return (baseHash + rotation) % ((leafSize + 1) / 2);
            }
        }
};

template <size_t leafSize>
class CandidateList {
    private:
        static constexpr uint64_t ALL_SET = (leafSize >= 127) ? ~0ul : (1ul << ((leafSize + 1) / 2)) - 1;
        std::vector<uint64_t> candidates;
    public:
        explicit CandidateList(size_t expectedNumSeeds) {
            candidates.reserve(expectedNumSeeds);
        }

        inline void add(size_t seed, uint64_t mask) {
            assert(seed == candidates.size());
            candidates.emplace_back(mask);
        }

        struct IteratorType {
            const CandidateList &candidateList;
            size_t currentIdx = -1;
            size_t size = 0;
            uint64_t filterMask;
            bool isEnd;

            IteratorType(const CandidateList &candidateList, uint64_t filterMask, bool isEnd)
                : candidateList(candidateList), filterMask(filterMask), isEnd(isEnd) {
                size = candidateList.candidates.size();
                if (!isEnd) {
                    operator++(); // Initialized with -1, go to first actual item
                }
            }

            inline bool operator!=(IteratorType rhs) const {
                return isEnd != rhs.isEnd;
            }

            inline std::pair<size_t, uint64_t> operator*() const {
                uint64_t mask = candidateList.candidates.at(currentIdx);
                return std::make_pair(currentIdx, mask);
            }

            inline IteratorType& operator++() {
                ++currentIdx;
                while (currentIdx < size) {
                    if ((candidateList.candidates[currentIdx] | filterMask) == ALL_SET) {
                        return *this;
                    }
                    ++currentIdx;
                }
                isEnd = true;
                return *this;
            }
        };

        struct FilteredListType {
            const CandidateList &candidateList;
            uint64_t mask;

            FilteredListType(CandidateList &candidateList, uint64_t mask)
                : candidateList(candidateList), mask(mask) {
            }

            inline IteratorType begin() const {
                return IteratorType(candidateList, mask, false);
            }

            inline IteratorType end() const {
                return IteratorType(candidateList, mask, true);
            }
        };

        FilteredListType filter(uint64_t mask) {
            return FilteredListType(*this, mask);
        }
};

template <size_t leafSize>
class CandidateTree {
    private:
        static constexpr uint64_t ALL_SET = (leafSize >= 127) ? ~0ul : (1ul << ((leafSize + 1) / 2)) - 1;
        static constexpr uint64_t BUCKET_MASK = 0b11111;
        std::array<std::vector<uint64_t>, BUCKET_MASK + 1> candidateMasks;
        std::array<std::vector<uint64_t>, BUCKET_MASK + 1> candidateSeeds;
    public:
        explicit CandidateTree(size_t expectedNumSeeds) {
            size_t toReserve = expectedNumSeeds / candidateSeeds.size();
            for (size_t i = 0; i < candidateSeeds.size(); i++) {
                candidateSeeds[i].reserve(toReserve);
                candidateMasks[i].reserve(toReserve);
            }
        }

        inline void add(size_t seed, uint64_t mask) {
            candidateMasks[mask & BUCKET_MASK].emplace_back(mask);
            candidateSeeds[mask & BUCKET_MASK].emplace_back(seed);
        }

        struct IteratorType {
            const CandidateTree &candidateTree;
            size_t currentBucket = 0;
            size_t currentIdx = -1;
            uint64_t filterMask;
            uint64_t filterMaskRestrictedToBucketMask;
            bool isEnd;
            size_t currentBucketSize;

            inline void nextRelevantBucket() {
                while ((filterMaskRestrictedToBucketMask | currentBucket) < BUCKET_MASK) {
                    currentBucket++;
                }
                if ((filterMaskRestrictedToBucketMask | currentBucket) < BUCKET_MASK && currentBucket <= BUCKET_MASK) {
                    currentBucket++;
                }
            }

            IteratorType(const CandidateTree &candidateTree, uint64_t filterMask, bool isEnd)
                    : candidateTree(candidateTree), filterMask(filterMask),
                      filterMaskRestrictedToBucketMask(filterMask & BUCKET_MASK), isEnd(isEnd) {
                if (!isEnd) {
                    nextRelevantBucket();
                    currentBucketSize = candidateTree.candidateMasks[currentBucket].size();
                    operator++(); // Initialized with -1, go to first actual item
                }
            }

            inline bool operator!=(IteratorType rhs) const {
                return isEnd != rhs.isEnd;
            }

            inline std::pair<size_t, uint64_t> operator*() const {
                return std::make_pair(candidateTree.candidateSeeds[currentBucket][currentIdx],
                                      candidateTree.candidateMasks[currentBucket][currentIdx]);
            }

            inline IteratorType& operator++() {
                ++currentIdx;
                while (currentBucket <= BUCKET_MASK) {
                    while (currentIdx < currentBucketSize) {
                        if ((candidateTree.candidateMasks[currentBucket][currentIdx] | filterMask) == ALL_SET) {
                            return *this;
                        }
                        ++currentIdx;
                    }
                    currentBucket++;
                    nextRelevantBucket();
                    currentBucketSize = candidateTree.candidateMasks[currentBucket].size();
                    currentIdx = 0;
                }
                isEnd = true;
                return *this;
            }
        };

        struct FilteredListType {
            const CandidateTree &candidateTree;
            uint64_t mask;

            FilteredListType(CandidateTree &candidateTree, uint64_t mask)
                    : candidateTree(candidateTree), mask(mask) {
            }

            inline IteratorType begin() const {
                return IteratorType(candidateTree, mask, false);
            }

            inline IteratorType end() const {
                return IteratorType(candidateTree, mask, true);
            }
        };

        FilteredListType filter(uint64_t mask) {
            return FilteredListType(*this, mask);
        }
};

template <template<size_t> typename CandidateList, size_t leafSize, bool filter = false>
class QuadSplitCandidateFinder {
    private:
        static constexpr double E_HALF = 1.359140914;
        std::array<uint64_t, leafSize> keys = { 0 };
        uint64_t sizeSetA = 0;
        size_t currentSeed = -1;
        CandidateList<leafSize> candidatesA;
        CandidateList<leafSize> candidatesB;
        std::vector<SeedCache<leafSize, filter>> extractedCandidates;
    public:
        explicit QuadSplitCandidateFinder(const std::vector<uint64_t> &keysIn)
                : candidatesA(std::pow(E_HALF, keysIn.size() / 2)), candidatesB(std::pow(E_HALF, keysIn.size() / 2)) {
            for (size_t i = 0; i < leafSize; i++) {
                if ((keysIn[i] & 1) == 0) {
                    keys.at(sizeSetA) = keysIn[i];
                    sizeSetA++;
                } else {
                    keys.at(leafSize - i + sizeSetA - 1) = keysIn[i];
                }
            }
        }

        inline SeedCache<leafSize, filter> next() {
            while (extractedCandidates.empty()) {
                prepareNextSeed();
            }
            SeedCache<leafSize, filter> seedCache = extractedCandidates.back(); // Smallest seed is in the back
            extractedCandidates.pop_back();
            return seedCache;
        }

        inline SeedCache<leafSize, filter> makeCache(size_t seedA, size_t seedB) {
            SeedCache<leafSize, filter> cache = {};
            cache.seed = pairElegant(seedA, seedB);
            for (size_t i = 0; i < sizeSetA; i++) {
                uint64_t hash = util::remix(keys[i] + seedA);
                cache.hashes[i] = util::fastrange64(hash, (leafSize + 1) / 2);
            }
            for (size_t i = sizeSetA; i < leafSize; i++) {
                uint64_t hash = util::remix(keys[i] + seedB);
                cache.hashes[i] = util::fastrange64(hash, (leafSize + 1) / 2);
            }
            if constexpr (filter) {
                calculateIsolatedVertices(cache);
            }
            return cache;
        }

        void prepareNextSeed() {
            currentSeed++;

            uint64_t takenA = 0;
            uint64_t takenB = 0;
            for (size_t i = 0; i < sizeSetA; i++) {
                uint64_t hash = util::remix(keys[i] + currentSeed);
                takenA |= 1ul << util::fastrange64(hash, (leafSize + 1) / 2);
            }
            for (size_t i = sizeSetA; i < leafSize; i++) {
                uint64_t hash = util::remix(keys[i] + currentSeed);
                takenB |= 1ul << util::fastrange64(hash, (leafSize + 1) / 2);
            }

            for (const auto [candidateSeed, candidateMask] : candidatesA.filter(takenB)) {
                if (isFloatAccurateElegant(candidateSeed, currentSeed)) {
                    extractedCandidates.push_back(makeCache(candidateSeed, currentSeed));
                } else {
                    std::cout << "Skipped seed because of floating point inaccuracy" << std::endl;
                }
            }
            candidatesA.add(currentSeed, takenA); // add after iterating, so we don't test the same seed combination twice
            candidatesB.add(currentSeed, takenB);
            for (const auto [candidateSeed, candidateMask] : candidatesB.filter(takenA)) {
                if (isFloatAccurateElegant(currentSeed, candidateSeed)) {
                    extractedCandidates.push_back(makeCache(currentSeed, candidateSeed));
                } else {
                    std::cout << "Skipped seed because of floating point inaccuracy" << std::endl;
                }
            }
            if (extractedCandidates.size() > 1) {
                // Smallest seed is in the back
                std::sort(extractedCandidates.begin(), extractedCandidates.end(),
                          [](const SeedCache<leafSize, filter> &a, const SeedCache<leafSize, filter> &b) {
                              return a.seed > b.seed; });
            }
        }

        static size_t hash(uint64_t key, uint64_t seed) {
            auto [seedA, seedB] = unpairElegant(seed);
            if ((key & 1) == 0) {
                return util::fastrange64(util::remix(key + seedA), (leafSize + 1) / 2);
            } else {
                return util::fastrange64(util::remix(key + seedB), (leafSize + 1) / 2);
            }
        }
};

template <size_t leafSize, bool filter>
using QuadSplitCandidateFinderList = QuadSplitCandidateFinder<CandidateList, leafSize, filter>;

template <size_t leafSize, bool filter>
using QuadSplitCandidateFinderTree = QuadSplitCandidateFinder<CandidateTree, leafSize, filter>;

/**
 * ShockHash2 base case.
 * Note that while this can be used with uneven leaf sizes, it achieves suboptimal space and time.
 */
template <size_t leafSize, bool filter = false,
        template<size_t, bool> typename SeedCandidateFinder = BasicSeedCandidateFinder>
class BijectionsShockHash2 {
    public:
        static inline size_t findSeed(const std::vector<uint64_t> &keys) {
            assert(keys.size() == leafSize);
            SeedCandidateFinder<leafSize, filter> seedCandidateFinder(keys);
            std::vector<SeedCache<leafSize, filter>> seedsCandidates;
            seedsCandidates.push_back(seedCandidateFinder.next());
            CuckooUnionFind unionFind(leafSize);
            while (true) {
                const SeedCache<leafSize, filter> newCandidate = seedCandidateFinder.next();
                SeedCache<leafSize, filter> newCandidateShifted = newCandidate;
                for (size_t i = 0; i < leafSize; i++) {
                    newCandidateShifted.hashes[i] += leafSize/2;
                }
                for (const SeedCache<leafSize, filter> &other : seedsCandidates) {
                    if constexpr (filter) {
                        if ((other.isolatedVertices & newCandidate.isolatedVertices) != 0) {
                            continue;
                        }
                    }
                    size_t i = 0;
                    unionFind.clear();
                    for (; i < leafSize; i++) {
                        size_t end1 = other.hashes[i];
                        size_t end2 = newCandidateShifted.hashes[i];
                        if (!unionFind.unionIsStillPseudoforest(end1, end2)) {
                            break;
                        }
                    }
                    if (i != leafSize) {
                        continue;
                    }

                    // Found working seed!
                    uint64_t seed1 = newCandidate.seed;
                    uint64_t seed2 = other.seed;
                    uint64_t fullSeed = pairTriangular(seed1, seed2);

                    if (!isFloatAccurateTriangular(seed1, seed2)) {
                        std::cout << "Skipped seed because of floating point inaccuracy" << std::endl;
                        continue;
                    }
                    return fullSeed;
                }
                seedsCandidates.push_back(newCandidate);
            }
            return 0;
        }

        static inline void constructRetrieval(const std::vector<uint64_t> &keys, size_t seed,
                                       std::vector<std::pair<uint64_t, uint8_t>> &retrieval) {
            auto [seed1, seed2] = unpairTriangular(seed);
            shockhash::TinyBinaryCuckooHashTable table(leafSize);
            for (uint64_t key : keys) {
                table.prepare(shockhash::HashedKey(key));
            }
            table.clearPlacement();
            for (size_t k = 0; k < leafSize; k++) {
                shockhash::TinyBinaryCuckooHashTable::CandidateCells candidateCells;
                candidateCells.cell1 = SeedCandidateFinder<leafSize, filter>::hash(table.heap[k].hash.mhc, seed1) + leafSize/2;
                candidateCells.cell2 = SeedCandidateFinder<leafSize, filter>::hash(table.heap[k].hash.mhc, seed2);
                if (!table.insert(&table.heap[k], candidateCells)) {
                    throw std::logic_error("Should be possible to construct");
                }
            }
            for (size_t k = 0; k < leafSize; k++) {
                size_t candidate2 = SeedCandidateFinder<leafSize, filter>::hash(table.heap[k].hash.mhc, seed2);
                if (table.cells[candidate2] == &table.heap[k]) {
                    retrieval.emplace_back(table.heap[k].hash.mhc, 1);
                } else {
                    retrieval.emplace_back(table.heap[k].hash.mhc, 0);
                }
            }
        }

        static inline size_t hash(size_t seed, uint64_t key, uint64_t retrieved) {
            auto [seed1, seed2] = unpairTriangular(seed);
            size_t result;
            if (retrieved == 0) {
                result = SeedCandidateFinder<leafSize, filter>::hash(key, seed1) + leafSize/2;
            } else {
                result = SeedCandidateFinder<leafSize, filter>::hash(key, seed2);
            }
            assert(result <= leafSize);
            return result;
        }

        static void verify(size_t seed, const std::vector<uint64_t> &keys,
                    std::vector<std::pair<uint64_t, uint8_t>> &retrieval) {
            std::vector<bool> taken(leafSize, false);
            for (uint64_t key : keys) {
                size_t retrieved = ~0u;
                for (auto &retr : retrieval) {
                    if (retr.first == key) {
                        retrieved = retr.second;
                    }
                }
                if (retrieved == ~0u) {
                    throw std::logic_error("Not in retrieval");
                }
                size_t hashValue = hash(seed, key, retrieved);
                if (taken[hashValue]) {
                    throw std::logic_error("Collision");
                }
                taken[hashValue] = true;
            }
        }

        inline double calculateBijection(const std::vector<uint64_t> &keys) {
            size_t seed = findSeed(keys);
            // Begin: Validity check
            #ifndef NDEBUG
                std::vector<std::pair<uint64_t, uint8_t>> retrieval;
                constructRetrieval(keys, seed, retrieval);
                verify(seed, keys, retrieval);
            #endif
            // End: Validity check
            return ceil(log2(seed + 1));
        }
};
