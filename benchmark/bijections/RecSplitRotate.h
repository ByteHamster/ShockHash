#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <sux/support/common.hpp>
#include <bytehamster/util/MurmurHash64.h>
#include <bytehamster/util/Function.h>

/**
 * Alternative method of finding bijections that rotates two halves.
 */
template <size_t leafSize, bool useLookupTable>
class BijectionsRotate {
    public:
        // Only used when useLookupTable is true.
        // normalize_rotations[x] contains i such that rotate(x, i) is minimal.
        // For two numbers (x,y) that can be rotated to match, we get:
        // rotate(x, normalize_rotations[x]) == rotate(y, normalize_rotations[y]).
        // To calculate how many times to rotate y to match x, we can subtract the number of rotations.
        std::vector<uint8_t> normalize_rotations;

        static inline constexpr uint32_t rotate(uint32_t val, uint8_t x) {
            return ((val << x) | (val >> (leafSize - x))) & ((1u<<leafSize) - 1);
        }

        // Temporary variables
        std::vector<uint64_t> itemsLeft;
        std::vector<uint64_t> itemsRight;

        BijectionsRotate() {
            itemsLeft.resize(leafSize);
            itemsRight.resize(leafSize);
            if constexpr (useLookupTable) {
                // This can probably be more efficient and calculated at compile time
                normalize_rotations.resize((1ul << leafSize) + 1);
                for (uint32_t pos = 1; pos != uint32_t(1ul << leafSize); pos++) {
                    uint32_t x = pos;
                    uint32_t min = x;
                    normalize_rotations[pos] = 0;
                    for (size_t i = 0; i < leafSize; i++) {
                        if (x < min) {
                            min = x;
                            normalize_rotations[pos] = i;
                        }
                        x = rotate(x, 1);
                    }
                }
            }
        }

        inline size_t calculateBijection(std::vector<uint64_t> &keys) {
            static_assert(leafSize <= 32, "Using uint32_t, so 32 is the maximum leaf size");

            // Split objects into two groups ("left" and "right")
            size_t itemsLeftCount = 0;
            size_t itemsRightCount = 0;
            for (int i = 0; i < leafSize; i++) {
                bool isLeft = (bytehamster::util::remix(keys[i] - 1) % 2) == 0;
                if (isLeft) {
                    itemsLeft[itemsLeftCount] = keys[i];
                    itemsLeftCount++;
                } else {
                    itemsRight[itemsRightCount] = keys[i];
                    itemsRightCount++;
                }
            }

            constexpr uint32_t found = uint32_t(1<<leafSize) - 1;
            for (size_t x = 0; true; x += leafSize) {
                uint32_t maskLeft = 0;
                uint32_t maskRight = 0;
                for (size_t i = 0; i < itemsLeftCount; i++) {
                    maskLeft |= uint32_t(1) << bytehamster::util::fastrange16(bytehamster::util::remix(itemsLeft[i] + x), leafSize);
                }
                if ((size_t) sux::nu(maskLeft) != itemsLeftCount) {
                    continue; // Collisions in left part
                }
                for (size_t i1 = 0; i1 < itemsRightCount; i1++) {
                    maskRight |= uint32_t(1) << bytehamster::util::fastrange16(bytehamster::util::remix(itemsRight[i1] + x), leafSize);
                }
                if ((size_t) sux::nu(maskRight) != itemsRightCount) {
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
