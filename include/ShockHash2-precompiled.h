#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace shockhash {
    size_t shockhash2query(size_t param, size_t seed, uint64_t key, size_t retrieved);
    size_t shockhash2construct(size_t param, std::vector<uint64_t> &leafKeys,
                               std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput);
}