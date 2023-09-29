#pragma once

#include <cmath>
#include <TinyBinaryCuckooHashTable.h>

/**
 * ShockHash base case
 */
template <size_t leafSize>
class BijectionsShockHash1 {
        static_assert(leafSize <= 64);
    public:
        inline double calculateBijection(std::vector<uint64_t> &keys) {
            shockhash::TinyBinaryCuckooHashTable tinyBinaryCuckooHashTable(leafSize);
            for (size_t i = 0; i < leafSize; i++) {
                tinyBinaryCuckooHashTable.prepare(shockhash::HashedKey(keys[i]));
            }
            constexpr uint64_t allSet = (leafSize == 64) ? (~0ul) : (1ul << leafSize) - 1;
            uint64_t mask = 0;
            size_t x = 0;
            for (;;) {
                for (;;) {
                    mask = 0;
                    for (size_t i = 0; i < leafSize; i++) {
                        auto hash = shockhash::TinyBinaryCuckooHashTable::getCandidateCells(shockhash::HashedKey(keys[i]), x, leafSize);
                        mask |= (1ul << hash.cell1);
                        mask |= (1ul << hash.cell2);
                    }
                    if (mask == allSet) break;
                    x++;
                }
                if (tinyBinaryCuckooHashTable.construct(x)) break;
                x++;
            }
            return ceil(log2(x + 1));
        }
};
