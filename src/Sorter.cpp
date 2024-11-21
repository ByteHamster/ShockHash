#include <ips2ra.hpp>
#include <sux/function/RecSplit.hpp>
#include "shockhash/Sorter.hpp"
#include "shockhash/ShockHash2FlatBase.h"

void shockhash::sort_hash128_t(sux::function::hash128_t *objects, size_t n, size_t threads) {
    ips2ra::parallel::sort(objects, objects + n, [&](const sux::function::hash128_t &a) { return a.first; }, threads);
}

void shockhash::sort_keyInfo(std::vector<shockhash::KeyInfo> &hashes) {
    ips2ra::sort(hashes.begin(), hashes.end(), [] (const KeyInfo &t) { return uint64_t(t.bucket) << 32 | t.threshold; });
}
