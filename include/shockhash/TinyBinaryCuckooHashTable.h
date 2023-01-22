#pragma once
#include <vector>
#include <cassert>
#include <queue>
#include "Function.h"
#include "MurmurHash64.h"
#include <cstring>

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
        uint8_t *unusedCells;
        size_t N;
        size_t M;
    private:
        size_t seed = 0;
        size_t numEntries = 0;
    public:
        explicit TinyBinaryCuckooHashTable(size_t N, size_t M) : N(N), M(M) {
            heap = new TableEntry[N];
            cells = new TableEntry*[M];
            unusedCells = new uint8_t[M];
        }

        ~TinyBinaryCuckooHashTable() {
            delete[] heap;
            delete[] cells;
            delete[] unusedCells;
        }

        void prepare(HashedKey hash) {
            assert(numEntries < N);
            heap[numEntries].hash = hash;
            numEntries++;
        }

        bool construct(size_t seed_) {
            seed = seed_;
            memset(unusedCells, 1, M);
            size_t differentHashes = 0;
            for (size_t i = 0; i < numEntries; i++) {
                Union64 hash = getCandidateCells(heap[i]);
                differentHashes += unusedCells[hash.halves.low];
                differentHashes += unusedCells[hash.halves.high];
                unusedCells[hash.halves.low] = 0;
                unusedCells[hash.halves.high] = 0;
            }
            if (differentHashes < N) {
                return false;
            }
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
            Union64 hash;
            hash.full = util::remix(key.mhc + seed);
            if (hashFunctionIndex == 0) {
                return util::fastrange32(hash.halves.high, range);
            } else {
                return util::fastrange32(hash.halves.low, range);
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

        inline Union64 getCandidateCells(TableEntry &entry) {
            Union64 hash;
            hash.full = util::remix(entry.hash.mhc + seed);
            hash.halves.high = util::fastrange32(hash.halves.high, M);
            hash.halves.low = util::fastrange32(hash.halves.low, M);
            return hash;
        }

        bool insert(TableEntry *entry) {
            Union64 candidates = getCandidateCells(*entry);
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
} // Namespace sichash
