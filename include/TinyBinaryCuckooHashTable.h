#pragma once
#include <vector>
#include <cassert>
#include <queue>
#include "Function.h"
#include "MurmurHash64.h"
#include <cstring>
#include "UnionFind.h"
#ifdef SIMD
#include "SimdUtils.h"
#endif

namespace shockhash {
struct HashedKey {
    uint64_t mhc;

    HashedKey() {
        this->mhc = 0;
    }

    explicit HashedKey(uint64_t mhc) : mhc(mhc) {
    }

    explicit HashedKey(const std::string &element, uint32_t seed = 0) {
        uint64_t stringHash = util::MurmurHash64(element.data(), element.length());
        uint64_t modified = stringHash + seed;
        mhc = util::MurmurHash64(&modified, sizeof(uint64_t));
    }

    [[nodiscard]] inline uint64_t hash(int hashFunctionIndex, size_t range) const {
        return util::fastrange64(util::remix(mhc + hashFunctionIndex), range);
    }
};

/**
 * Tiny binary cuckoo hash table. Construction needs multiple tries before succeeding.
 */
class TinyBinaryCuckooHashTable {
    public:
        struct TableEntry {
            HashedKey hash;
            uint32_t candidateCellsXor = 0;
        };
        TableEntry *heap;
        TableEntry** cells;
        size_t N;
        size_t M;
        UnionFind unionFind;
    private:
        size_t seed = 0;
        size_t numEntries = 0;
    public:
        explicit TinyBinaryCuckooHashTable(size_t N, size_t M) : N(N), M(M), unionFind(M) {
            heap = new TableEntry[N];
            cells = new TableEntry*[M];
        }

        ~TinyBinaryCuckooHashTable() {
            delete[] heap;
            delete[] cells;
        }

        void prepare(HashedKey hash) {
            assert(numEntries < N);
            heap[numEntries].hash = hash;
            numEntries++;
        }

        bool construct(size_t seed_) {
            seed = seed_;

            // Filter based on unused cells
            assert(numEntries <= 64);
            /*uint64_t used = 0;
            for (size_t i = 0; i < numEntries; i++) {
                Union64 hash = getCandidateCells(heap[i].hash, seed, M);
                used |= (1ul << hash.halves.low) | (1ul << hash.halves.high);
            }
            if (std::popcount(used) != numEntries) {
                return false;
            }

            // Filter based on tree/pseudotree
            unionFind.clear();
            for (size_t i = 0; i < numEntries; i++) {
                Union64 hash = getCandidateCells(heap[i].hash, seed, M);
                if (!unionFind.unionIsStillPseudoforrest(hash.halves.high, hash.halves.low)) {
                    return false;
                }
            }*/

            // Actual cuckoo hashing, we know that this will succeed at this point
            memset(cells, 0, M * sizeof(void*)); // Fill with nullpointers
            for (size_t i = 0; i < numEntries; i++) {
                if (!insert(&heap[i])) {
                    return false;
                }
            }
            return true;
        }

        [[nodiscard]] size_t size() const {
            return numEntries;
        }

        static inline size_t hashToCell(HashedKey key, size_t seed, size_t range, size_t hashFunctionIndex) {
            CandidateCells hash = getCandidateCells(key, seed, range);
            if (hashFunctionIndex == 0) {
                return hash.cell2;
            } else {
                return hash.cell1;
            }
        }

        typedef struct {
            uint32_t cell1;
            uint32_t cell2;
        } CandidateCells;

        static inline CandidateCells getCandidateCells(const HashedKey key, size_t seed, size_t range) {
            uint64_t remixed = util::remix(key.mhc + seed);
            const uint32_t hash1 = util::fastrange32(remixed, range / 2);
            const uint32_t hash2 = util::fastrange32(remixed >> 32, (range + 1) / 2) + range / 2;
            return {hash1, hash2};
        }

#ifdef SIMD
        typedef struct {
            Vec4x64ui cell1;
            Vec4x64ui cell2;
        } CandidateCellsSIMD;

        static inline CandidateCellsSIMD getCandidateCellsSIMD(Vec4x64ui valueAndSeed, uint64_t range) {
            const Vec4x64ui remixed = remixV(valueAndSeed);
            const Vec4x64ui hash1 = remap32V(remixed, range / 2);
            const Vec4x64ui hash2 = remap32V(remixed >> 32, (range + 1) / 2) + range / 2;

            assert(TinyBinaryCuckooHashTable::getCandidateCells(HashedKey(valueAndSeed.extract(0)), 0, range).cell1 == hash1.extract(0));
            assert(TinyBinaryCuckooHashTable::getCandidateCells(HashedKey(valueAndSeed.extract(0)), 0, range).cell2 == hash2.extract(0));

            return {hash1, hash2};
        }
#endif

    private:

        bool insert(TableEntry *entry) {
            CandidateCells candidates = getCandidateCells(entry->hash, seed, M);
            entry->candidateCellsXor = candidates.cell1 ^ candidates.cell2;
            if (cells[candidates.cell1] == nullptr) {
                cells[candidates.cell1] = entry;
                return true;
            }
            if (cells[candidates.cell2] == nullptr) {
                cells[candidates.cell2] = entry;
                return true;
            }
            uint32_t currentCell = candidates.cell1;

            size_t tries = 0;
            while (tries < M) {
                uint32_t alternativeCell = entry->candidateCellsXor ^ currentCell;
                std::swap(entry, cells[alternativeCell]);
                if (entry == nullptr) {
                    return true;
                }
                currentCell = alternativeCell;
                tries++;
            }
            return false;
        }
};
} // Namespace shockhash
