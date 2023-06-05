#include <sux/function/RecSplit.hpp>

namespace shockhash {
    /**
     * Moving sorting to its own compile unit brings compile times from minutes to seconds.
     */
    void sort_hash128_t(sux::function::hash128_t *objects, size_t n);
}
