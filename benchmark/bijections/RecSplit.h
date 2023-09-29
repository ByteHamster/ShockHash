#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <sux/support/common.hpp>
#include <MurmurHash64.h>
#include <Function.h>

static const int MAX_LEAF_SIZE = 24;
#define DO_NOT_OPTIMIZE(value) asm volatile ("" : : "r,m"(value) : "memory")

/**
 * Copy of RecSplit's bijection code.
 */
template <size_t leafSize>
class BijectionsRecSplit {
    public:
        static constexpr std::array<uint8_t, MAX_LEAF_SIZE> fill_bij_midstop() {
            std::array<uint8_t, MAX_LEAF_SIZE> memo{0};
            for (int s = 0; s < MAX_LEAF_SIZE; ++s) memo[s] = s < (int)ceil(2 * sqrt(s)) ? s : (int)ceil(2 * sqrt(s));
            return memo;
        }
        static constexpr std::array<uint8_t, MAX_LEAF_SIZE> bij_midstop = fill_bij_midstop();

        inline double calculateBijection(std::vector<uint64_t> &keys) {
            uint32_t mask;
            int x = 0;
            const uint32_t found = (1 << leafSize) - 1;

            if constexpr (leafSize <= 8) {
                for (;;) {
                    mask = 0;
                    for (size_t i = 0; i < leafSize; i++) {
                        mask |= uint32_t(1) << util::fastrange16(util::remix(keys[i] + x), leafSize);
                    }
                    if (mask == found) {
                        return ceil(log2(x + 1));
                    }
                    x++;
                }
            } else {
                const size_t midstop = bij_midstop[leafSize];
                for (;;) {
                    mask = 0;
                    size_t i;
                    for (i = 0; i < midstop; i++) {
                        mask |= uint32_t(1) << util::fastrange16(util::remix(keys[i] + x), leafSize);
                    }
                    if (sux::nu(mask) == midstop) {
                        for (; i < leafSize; i++) {
                            mask |= uint32_t(1) << util::fastrange16(util::remix(keys[i] + x), leafSize);
                        }
                        if (mask == found) {
                            return ceil(log2(x + 1));
                        }
                    }
                    x++;
                }
            }
        }
};
