#pragma once
#include <vector>
#include <cassert>
#include <queue>
#include "Function.h"
#include "MurmurHash64.h"
#include <cstring>
#include "UnionFind.h"

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
            uint64_t used = 0;
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
            }

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
            Union64 hash = getCandidateCells(key, seed, range);
            if (hashFunctionIndex == 0) {
                return hash.halves.high;
            } else {
                return hash.halves.low;
            }
        }

    private:
        typedef union {
            struct {
                uint32_t low;
                uint32_t high;
            } halves;
            uint64_t full;
        } Union64;

        static inline Union64 getCandidateCells(HashedKey &key, size_t seed, size_t range) {
            Union64 hash;
            hash.full = util::remix(key.mhc + seed);
            hash.halves.high = util::fastrange32(hash.halves.high, range/2);
            hash.halves.low = util::fastrange32(hash.halves.low, (range+1)/2) + range/2;
            return hash;
        }

        bool insert(TableEntry *entry) {
            Union64 candidates = getCandidateCells(entry->hash, seed, M);
            entry->candidateCellsXor = candidates.halves.high ^ candidates.halves.low;
            if (cells[candidates.halves.high] == nullptr) {
                cells[candidates.halves.high] = entry;
                return true;
            }
            if (cells[candidates.halves.low] == nullptr) {
                cells[candidates.halves.low] = entry;
                return true;
            }
            uint32_t currentCell = candidates.halves.low;

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
