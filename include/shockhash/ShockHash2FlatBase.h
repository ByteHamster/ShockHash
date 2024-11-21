#pragma once
#include <cstdint>

namespace shockhash {
struct KeyInfo {
    uint64_t mhc;
    uint32_t bucket;
    uint32_t threshold;
};
}
