#include <Sorter.hpp>
#include <ips2ra.hpp>
#include <sux/function/RecSplit.hpp>

void shockhash::sort_hash128_t(sux::function::hash128_t *objects, size_t n) {
    ips2ra::sort(objects, objects + n, [&](const sux::function::hash128_t &a) { return a.first; });
}
