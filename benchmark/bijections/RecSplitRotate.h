#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <sux/support/common.hpp>
#include <MurmurHash64.h>
#include <Function.h>

/**
 * Alternative method of finding bijections that rotates two halves.
 */
template <uint8_t leafSize, bool useLookupTable>
class BijectionsRotate {
    public:
        // Only used when useLookupTable is true.
        // normalize_rotations[x] contains i such that rotate(x, i) is minimal.
        // For two numbers (x,y) that can be rotated to match, we get:
        // rotate(x, normalize_rotations[x]) == rotate(y, normalize_rotations[y]).
        // To calculate how many times to rotate y to match x, we can subtract the number of rotations.
        uint8_t normalize_rotations[1<<leafSize] = {0};

        static inline constexpr uint16_t rotate(uint16_t val, uint8_t x) {
            return ((val << x) | (val >> (leafSize - x))) & ((1<<leafSize) - 1);
        }

        // Temporary variables
        uint64_t itemsLeft[leafSize] = {0};
        uint64_t itemsRight[leafSize] = {0};

        BijectionsRotate() {
            if constexpr (useLookupTable) {
                // This can probably be more efficient and calculated at compile time
                for (uint16_t pos = 1; pos != uint16_t(1<<leafSize); pos++) {
                    uint16_t x = pos;
                    uint16_t min = x;
                    normalize_rotations[pos] = 0;
                    for (uint8_t i = 0; i < leafSize; i++) {
                        x = rotate(x, 1);
                        if (x < min) {
                            min = x;
                            normalize_rotations[pos] = i + 1;
                        }
                    }
                }
            }
        }

        inline double calculateBijection(std::vector<uint64_t> &keys) {
            static_assert(leafSize <= 32, "Using uint32_t, so 32 is the maximum leaf size");

            // Split objects into two groups ("left" and "right")
            int itemsLeftCount = 0;
            int itemsRightCount = 0;
            for (int i = 0; i < leafSize; i++) {
                bool isLeft = (util::remix(keys[i] - 1) % 2) == 0;
                if (isLeft) {
                    itemsLeft[itemsLeftCount] = keys[i];
                    itemsLeftCount++;
                } else {
                    itemsRight[itemsRightCount] = keys[i];
                    itemsRightCount++;
                }
            }

            constexpr uint32_t found = uint32_t(1<<leafSize) - 1;
            for (int x = 0; true; x += leafSize) {
                uint32_t maskLeft = 0;
                uint32_t maskRight = 0;
                for (size_t i = 0; i < itemsLeftCount; i++) {
                    maskLeft |= uint32_t(1) << util::fastrange16(util::remix(itemsLeft[i] + x), leafSize);
                }
                if (sux::nu(maskLeft) != itemsLeftCount) {
                    continue; // Collisions in left part
                }
                for (size_t i1 = 0; i1 < itemsRightCount; i1++) {
                    maskRight |= uint32_t(1) << util::fastrange16(util::remix(itemsRight[i1] + x), leafSize);
                }
                if (sux::nu(maskRight) != itemsRightCount) {
                    continue; // Collisions in right part
                }
                // Try to rotate right part to see if both together form a bijection
                if constexpr (useLookupTable) {
                    uint8_t rotations = (normalize_rotations[uint32_t(~maskRight) & found]
                                         - normalize_rotations[maskLeft] + leafSize) % leafSize;
                    if (found == (maskLeft | rotate(maskRight, rotations))) {
                        return x + rotations;
                    }
                } else {
                    for (int rotations = 0; rotations < leafSize; rotations++) {
                        if (found == (maskLeft | maskRight)) {
                            return x + rotations;
                        }
                        maskRight = rotate(maskRight, 1);
                    }
                }
            }
        }

        static std::string name() {
            return std::string("RecSplitRotate") + (useLookupTable ? "Lookup" : "");
        }
};
