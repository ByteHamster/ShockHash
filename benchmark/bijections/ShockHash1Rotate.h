#pragma once

#include <cmath>
#include "shockhash/TinyBinaryCuckooHashTable.h"

/**
 * ShockHash base case
 */
template <size_t leafSize>
class BijectionsShockHash1Rotate {
        static_assert(leafSize <= 64);
    public:
        static constexpr uint64_t rotate(size_t l, uint64_t val, uint32_t x) {
            return ((val << x) | (val >> (l - x))) & ((1ul << l) - 1);
        }

        inline size_t calculateBijection(std::vector<uint64_t> &keysInput) {
            size_t x = 0;
            constexpr uint64_t allSet = (leafSize == 64) ? (~0ul) : (1ul << leafSize) - 1;
            size_t r = 0;
            size_t keysGroupA = 0;
            size_t indexB = leafSize - 1;
            shockhash::HashedKey keys[leafSize];
            shockhash::TinyBinaryCuckooHashTable::CandidateCells candidateCellsCache[leafSize];
            shockhash::TinyBinaryCuckooHashTable tinyBinaryCuckooHashTable(leafSize);
            for (size_t i = 0; i < leafSize; i++) {
                auto key = shockhash::HashedKey(keysInput[i]);
                if ((key.mhc & 1) == 0) {
                    keys[keysGroupA] = key;
                    keysGroupA++;
                } else {
                    keys[indexB] = key;
                    indexB--;
                }
            }
            for (size_t i = 0; i < leafSize; i++) {
                tinyBinaryCuckooHashTable.prepare(keys[i]);
            }
            uint64_t a = 0;
            uint64_t b = 0;
            for (;;x++) {
                a = 0;
                for (size_t i = 0; i < keysGroupA; i++) {
                    auto candidateCells = shockhash::TinyBinaryCuckooHashTable::getCandidateCells<leafSize>(keys[i], x);
                    candidateCellsCache[i] = candidateCells;
                    uint64_t candidatePowers = (1ull << candidateCells.cell1) | (1ull << candidateCells.cell2);
                    a |= candidatePowers;
                }
                b = 0;
                for (size_t i = keysGroupA; i < leafSize; i++) {
                    auto candidateCells = shockhash::TinyBinaryCuckooHashTable::getCandidateCells<leafSize>(keys[i], x);
                    candidateCellsCache[i] = candidateCells;
                    uint64_t candidatePowers = (1ull << candidateCells.cell1) | (1ull << candidateCells.cell2);
                    b |= candidatePowers;
                }
                for (r = 0; r < leafSize; r++) {
                    if ((a | rotate(leafSize, b, r)) != allSet) {
                        continue;
                    }
                    tinyBinaryCuckooHashTable.clearPlacement();
                    size_t i = 0;
                    for (i = 0; i < leafSize; i++) {
                        auto candidateCells = candidateCellsCache[i];
                        if ((tinyBinaryCuckooHashTable.heap[i].hash.mhc & 1) == 1) {
                            // Set B
                            candidateCells.cell1 = (candidateCells.cell1 + r) % leafSize;
                            candidateCells.cell2 = (candidateCells.cell2 + r) % leafSize;
                        }
                        if (!tinyBinaryCuckooHashTable.insert(&tinyBinaryCuckooHashTable.heap[i], candidateCells)) {
                            break;
                        }
                    }
                    if (i == leafSize) {
                        return x * leafSize + r;
                    }
                }
            }
        }

        static std::string name() {
            return "ShockHashRotate";
        }
};
