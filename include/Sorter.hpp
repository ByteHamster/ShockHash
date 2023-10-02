#pragma once

#include <sux/function/RecSplit.hpp>
#include <vector>

namespace shockhash {
    /**
     * Moving sorting to its own compile unit brings compile times from minutes to seconds.
     */
    void sort_hash128_t(sux::function::hash128_t *objects, size_t n, size_t threads);

    void parallelPartition(sux::function::hash128_t *input, std::vector<uint64_t> &sorted,
                           std::vector<uint64_t> &bucket_size_acc, size_t num_threads,
                           size_t keys_count, size_t nbuckets) {
        assert(sorted.size() >= keys_count);
        assert(bucket_size_acc.size() == nbuckets + 1);

        sort_hash128_t(input, keys_count, num_threads);

        // For most reasonable input sizes, doing this sequentially is faster, see GpuRecSplit
        size_t i = 0;
        const sux::function::hash128_t *it = input;
        const sux::function::hash128_t *end = input + keys_count;
        for (size_t bucket = 0; bucket < nbuckets; bucket++) {
            bucket_size_acc.at(bucket) = i;
            while (sux::remap128(it->first, nbuckets) == bucket && it != end) {
                sorted[i] = it->second;
                i++;
                it++;
            }
        }
        bucket_size_acc[nbuckets] = keys_count;
    }
}
