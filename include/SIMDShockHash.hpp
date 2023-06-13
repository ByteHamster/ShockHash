#pragma once
/*
 * Based on RecSplit, Copyright (C) 2019-2020 Emmanuel Esposito and Sebastiano Vigna
 * Enhanced to use overloaded cuckoo hash tables in the leaves.
 * For tiny space usages (~1.6 bit/object), ShockHash is faster than RecSplit.
 */

#include <ShockHash.h>
#include <sux/util/Vector.hpp>
#include <sux/function/DoubleEF.hpp>
#include <sux/function/RiceBitVector.hpp>
#include "SimdUtils.h"
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <utility>

#ifndef SIMD
#error Need to compile SIMDShockHash with -DSIMD
#endif

namespace shockhash {

using namespace std;
using namespace std::chrono;
using namespace sux::function;

#if defined(MORESTATS) && !defined(STATS)
#define STATS
#endif

#ifdef MORESTATS

#define MAX_LEVEL_TIME (20)

static constexpr double log2e = 1.44269504089;
static uint64_t num_bij_trials[MAX_LEAF_SIZE], num_split_trials;
static uint64_t num_bij_evals[MAX_LEAF_SIZE], num_split_evals;
static uint64_t bij_count[MAX_LEAF_SIZE], split_count;
static uint64_t expected_split_trials, expected_split_evals;
static uint64_t bij_unary, bij_fixed, bij_unary_golomb, bij_fixed_golomb;
static uint64_t split_unary, split_fixed, split_unary_golomb, split_fixed_golomb;
static uint64_t max_split_code, min_split_code, sum_split_codes;
static uint64_t max_bij_code, min_bij_code, sum_bij_codes;
static uint64_t sum_depths;
static uint64_t time_bij;
static uint64_t time_split[MAX_LEVEL_TIME];
size_t minsize = 0, maxsize = 0;
double ub_split_bits = 0, ub_bij_bits = 0;
double ub_split_evals = 0;
#endif

/**
  * 32-bit finalizer function in Austin Appleby's MurmurHash3 (https://github.com/aappleby/smhasher).
  */
FullVecUi inline remix32V(FullVecUi z) {
    z ^= z >> 16;
    z *= 0x85ebca6b;
    z ^= z >> 13;
    z *= 0xc2b2ae35;
    z ^= z >> 16;
    return z;
}

FullVecUi remap(FullVecUq x, FullVecUq y, uint32_t n) {
    constexpr int masklen = 16;
    constexpr uint32_t mask = (uint32_t(1) << masklen) - 1;
    //FullVecUi combined(compress(x >> 32, y >> 32) ^ compress(x, y)); // This is a bit slower than below
    const FullVecUi xx(x);
    const FullVecUi yy(y);
#ifdef SIMDRS_512_BIT
    FullVecUi combined = blend16<0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30>(xx, yy)
        ^ blend16<1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31>(xx, yy);
#else
    FullVecUi combined = blend8<0, 2, 4, 6, 8, 10, 12, 14>(xx, yy) ^ blend8<1, 3, 5, 7, 9, 11, 13, 15>(xx, yy);
#endif
    return ((remix32V(combined) & mask) * n) >> masklen;
}

/**
  * 32-bit finalizer function in Austin Appleby's MurmurHash3 (https://github.com/aappleby/smhasher).
  */
uint32_t remix32(uint32_t z) {
    z ^= z >> 16;
    z *= 0x85ebca6b;
    z ^= z >> 13;
    z *= 0xc2b2ae35;
    z ^= z >> 16;
    return z;
}

static FullVecUi powerOfTwo(FullVecUi x) {
#ifdef SIMDRS_512_BIT
    return _mm512_sllv_epi32(FullVecUi(1), x);
#else
    return _mm256_sllv_epi32(FullVecUi(1), x);
#endif
}

// The first 4 entries are 0, then there are the first 4 powers of 256
// and then a zero padding.
static constexpr array<uint32_t, 4 + MAX_FANOUT + FULL_VEC_32_COUNT> fill_aggr_level_count_lookup() {
    array<uint32_t, 4 + MAX_FANOUT + FULL_VEC_32_COUNT> memo{ 0 };
    for (size_t s = 0; s < 4; ++s) memo[s + 4] = 1 << (8 * s);
    return memo;
}

/**
 *
 * A class for storing minimal perfect hash functions. The template
 * parameter decides how large a leaf will be. Larger leaves imply
 * slower construction, but less space and faster evaluation.
 *
 * @tparam LEAF_SIZE the size of a leaf; typicals value range from 6 to 8
 * for fast, small maps, or up to 16 for very compact functions.
 * @tparam AT a type of memory allocation out of sux::util::AllocType.
 */

template <size_t LEAF_SIZE, sux::util::AllocType AT = sux::util::AllocType::MALLOC>
class SIMDShockHash {
    static_assert(LEAF_SIZE <= MAX_LEAF_SIZE);
    static constexpr size_t _leaf = LEAF_SIZE;
    static constexpr size_t lower_aggr = SplittingStrategy<LEAF_SIZE>::lower_aggr;
    static constexpr size_t upper_aggr = SplittingStrategy<LEAF_SIZE>::upper_aggr;

    // For each bucket size, the Golomb-Rice parameter (upper 8 bits) and the number of bits to
    // skip in the fixed part of the tree (lower 24 bits).
    static constexpr array<uint32_t, MAX_BUCKET_SIZE> memo = fill_golomb_rice<LEAF_SIZE>();
    static constexpr auto aggr_level_count_lookup = fill_aggr_level_count_lookup();

    size_t bucket_size;
    size_t nbuckets;
    size_t keys_count;
    RiceBitVector<AT> descriptors;
    DoubleEF<AT> ef;
    SimpleRibbon<1> *ribbon = nullptr;
    std::vector<std::pair<uint64_t, uint8_t>> ribbonInput;

  public:
    SIMDShockHash() {}

    /** Builds a SIMDShockHash instance using a given list of keys and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * @param keys a vector of strings.
     * @param bucket_size the desired bucket size; typical sizes go from
     * 100 to 2000, with smaller buckets giving slightly larger but faster
     * functions.
     */
    SIMDShockHash(const vector<string> &keys, const size_t bucket_size) {
        this->bucket_size = bucket_size;
        this->keys_count = keys.size();
        hash128_t *h = (hash128_t *)malloc(this->keys_count * sizeof(hash128_t));
        for (size_t i = 0; i < this->keys_count; ++i) {
            h[i] = first_hash(keys[i].c_str(), keys[i].size());
        }
        hash_gen(h);
        free(h);
    }

    /** Builds a SIMDShockHash instance using a given list of 128-bit hashes and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * Note that this constructor is mainly useful for benchmarking.
     * @param keys a vector of 128-bit hashes.
     * @param bucket_size the desired bucket size; typical sizes go from
     * 100 to 2000, with smaller buckets giving slightly larger but faster
     * functions.
     */
    SIMDShockHash(vector<hash128_t> &keys, const size_t bucket_size) {
        this->bucket_size = bucket_size;
        this->keys_count = keys.size();
        hash_gen(&keys[0]);
    }

    /** Builds a SIMDShockHash instance using a list of keys returned by a stream and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * @param input an open input stream returning a list of keys, one per line.
     * @param bucket_size the desired bucket size.
     */
    SIMDShockHash(ifstream& input, const size_t bucket_size) {
        this->bucket_size = bucket_size;
        vector<hash128_t> h;
        for(string key; getline(input, key);) h.push_back(first_hash(key.c_str(), key.size()));
        this->keys_count = h.size();
        hash_gen(&h[0]);
    }

  private:
    // Maps a 128-bit to a bucket using the first 64-bit half.
    inline uint64_t hash128_to_bucket(const hash128_t &hash) const { return remap128(hash.first, nbuckets); }

    uint64_t higherLevel(vector<uint64_t> &bucket, vector<uint64_t> &temp, size_t start, size_t m, size_t split, typename RiceBitVector<AT>::Builder &builder,
        vector<uint32_t> &unary, const int level) {
#ifdef MORESTATS
        auto start_time = high_resolution_clock::now();
#endif
        const auto end = start + m;
        assert(m > upper_aggr);
        assert(m <= MAX_BUCKET_SIZE);
        uint64_t x = start_seed[level];
#ifdef SIMDRS_512_BIT
        FullVecUq xVec(0, 1, 2, 3, 4, 5, 6, 7);
#else
        FullVecUq xVec(0, 1, 2, 3);
#endif
        xVec += x;

        FullVecUi counter;
        FullVecIb found_result;
        for (;;) {
            counter = 0;
            for (size_t i = start; i < end; i++) {
                const FullVecUq first = bucket[i] + xVec;
                // seems equal to "counter += remap(first, first + FULL_VEC_64_COUNT, m) < split" and faster than if_sub
                counter = if_add(remap(first, first + FULL_VEC_64_COUNT, m) < split, counter, uint32_t(-1));
#ifdef MORESTATS
                ++num_split_evals;
#endif
            }
            found_result = counter == -split; // -split since true is represented as -1 in vectors
            if (horizontal_or(found_result)) break;
            x += FULL_VEC_32_COUNT;
            xVec += FULL_VEC_32_COUNT;
        }

        const auto found_idx = horizontal_find_first(found_result);
        assert(found_idx != -1);
        x += found_idx;

        size_t count[2];
        count[0] = 0;
        count[1] = split;
        size_t i;
        for (i = start; i + FULL_VEC_32_COUNT <= end; i += FULL_VEC_32_COUNT) {
            FullVecUq first, second;
            first.load(&bucket[i]);
            second.load(&bucket[i + FULL_VEC_64_COUNT]);
            auto bits = to_bits(remap(first + x, second + x, m) >= split); // Fast for AVX-512
            for (size_t j = 0; j < FULL_VEC_32_COUNT; ++j)
                temp[count[(bits >> j) & 1]++] = bucket[i + j]; // TODO: Vectorize? Probably hard!
        }
        FullVecUq first, second;
        first.load(&bucket[i]);
        second.load(&bucket[i + FULL_VEC_64_COUNT]);
        auto bits = to_bits(remap(first + x, second + x, m) >= split);
        for (size_t j = 0; j < end - i; ++j)
            temp[count[(bits >> j) & 1]++] = bucket[i + j];
        copy(&temp[0], &(temp.data()[m]), &bucket[start]);

        x -= start_seed[level];
        const auto log2golomb = golomb_param(m);
        builder.appendFixed(x, log2golomb);
        unary.push_back(x >> log2golomb);

#ifdef MORESTATS
        time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
        return x;
    }

    template<uint32_t _MAX_FANOUT, uint32_t SPLIT, uint64_t SEED>
    uint64_t aggrLevel(vector<uint64_t> &bucket, vector<uint64_t> &temp, size_t start, size_t m, typename RiceBitVector<AT>::Builder &builder,
        vector<uint32_t> &unary, [[maybe_unused]] const int level) {
#ifdef MORESTATS
        auto start_time = high_resolution_clock::now();
#endif
        const auto end = start + m;
        uint64_t x = SEED;
        const uint32_t fanout = uint16_t(m + SPLIT - 1) / SPLIT;
        assert(m > LEAF_SIZE);
        assert(m <= _MAX_FANOUT * SPLIT);
        assert(fanout >= 2);
        assert(fanout <= _MAX_FANOUT);
        size_t i;
        static_assert(_MAX_FANOUT <= MAX_FANOUT, "_MAX_FANOUT must be at most MAX_FANOUT!");
        static_assert(SPLIT < 255, "SPLIT must be less than 255 for aggrLevel to work correctly!"
            "Note that less than 256 is not enough since an overflow may carry to the next count.");
#ifdef SIMDRS_512_BIT
        FullVecUq xVec(SEED, SEED + 1, SEED + 2, SEED + 3, SEED + 4, SEED + 5, SEED + 6, SEED + 7);
#else
        FullVecUq xVec(SEED, SEED + 1, SEED + 2, SEED + 3);
#endif
        FullVecIb found_result;
        if (fanout <= 5) {
            FullVecUi count;
            uint32_t found = SPLIT;
            for (i = 1; i < fanout - 1; ++i)
                found |= SPLIT << (8 * i);
            if (fanout <= 4)
                found |= (m - (fanout - 1) * SPLIT) << (8 * i);
            for (;;) {
                count = 0;
                for (size_t i = start; i < end; i++) {
                    const FullVecUq first = bucket[i] + xVec;
                    const FullVecUi remapped = remap(first, first + FULL_VEC_64_COUNT, m) / const_uint(SPLIT);
                    count += lookup<min(_MAX_FANOUT, 5)>(remapped, &aggr_level_count_lookup[4]);
#ifdef MORESTATS
                    ++num_split_evals;
#endif
                }
                found_result = count == found;
                if (horizontal_or(found_result)) break;
                x += FULL_VEC_32_COUNT;
                xVec += FULL_VEC_32_COUNT;
            }
        } else {
            FullVecUi count_low;
            FullVecUi count_high;
            constexpr uint32_t found_low = (SPLIT << 24) | (SPLIT << 16) | (SPLIT << 8) | SPLIT;
            uint32_t found_high = SPLIT;
            for (i = 1; i < fanout - 5; ++i)
                found_high |= SPLIT << (8 * i);
            if (fanout <= 8)
                found_high |= (m - (fanout - 1) * SPLIT) << (8 * i);
            for (;;) {
                count_low = count_high = 0;
                for (size_t i = start; i < end; i++) {
                    const FullVecUq first = bucket[i] + xVec;
                    const FullVecUi remapped = remap(first, first + FULL_VEC_64_COUNT, m) / const_uint(SPLIT);
                    count_low += lookup<_MAX_FANOUT>(remapped, &aggr_level_count_lookup[4]);
                    count_high += lookup<_MAX_FANOUT>(remapped, &aggr_level_count_lookup[0]);
#ifdef MORESTATS
                    ++num_split_evals;
#endif
                }
                found_result = count_low == found_low & count_high == found_high;
                if (horizontal_or(found_result)) break;
                x += FULL_VEC_32_COUNT;
                xVec += FULL_VEC_32_COUNT;
            }
        }

        const auto found_idx = horizontal_find_first(found_result);
        assert(found_idx != -1);
        x += found_idx;

        uint64_t *temp_c = &temp[m];
        for (size_t i = 0, c = 0; i < fanout; i++, c += SPLIT) temp_c[i] = c;
        uint32_t remapped[FULL_VEC_32_COUNT];
        for (i = start; i + FULL_VEC_32_COUNT <= end; i += FULL_VEC_32_COUNT) {
            FullVecUq first, second;
            first.load(&bucket[i]);
            second.load(&bucket[i + FULL_VEC_64_COUNT]);
            (remap(first + x, second + x, m) / const_uint(SPLIT)).store(remapped);
            for (size_t j = 0; j < FULL_VEC_32_COUNT; ++j)
                temp[temp_c[remapped[j]]++] = bucket[i + j]; // TODO: Vectorize? Probably hard!
        }
        FullVecUq first, second;
        first.load(&bucket[i]);
        second.load(&bucket[i + FULL_VEC_64_COUNT]);
        (remap(first + x, second + x, m) / const_uint(SPLIT)).store(remapped);
        for (size_t j = 0; j < end - i; ++j)
            temp[temp_c[remapped[j]]++] = bucket[i + j];
        copy(&temp[0], &(temp.data()[m]), &bucket[start]);

        x -= SEED;
        const auto log2golomb = golomb_param(m);
        builder.appendFixed(x, log2golomb);
        unary.push_back(x >> log2golomb);

#ifdef MORESTATS
        time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
        return x;
    }

    static FullVecUi rotateRight(FullVecUi val, uint32_t x) {
        return ((val << uint32_t(LEAF_SIZE - x)) | (val >> x)) & ((1 << LEAF_SIZE) - 1);
    }

    void leafLevel(vector<uint64_t> &bucket, size_t start, size_t m, typename RiceBitVector<AT>::Builder &builder,
            vector<uint32_t> &unary, [[maybe_unused]] const int level) {
        assert(m >= 2);
        assert(m <= LEAF_SIZE);
        const auto end = start + m;
        constexpr uint64_t SEED = start_seed[NUM_START_SEEDS - 1];
        uint64_t x = SEED;
#ifdef MORESTATS
        sum_depths += m * level;
        auto start_time = high_resolution_clock::now();
#endif
        FullVecUi mask;
        const uint32_t found = (1 << m) - 1;
#ifdef SIMDRS_512_BIT
        FullVecUq xVec(SEED, SEED + 1, SEED + 2, SEED + 3, SEED + 4, SEED + 5, SEED + 6, SEED + 7);
        constexpr size_t midstop_threshold = 12;
#else
        FullVecUq xVec(SEED, SEED + 1, SEED + 2, SEED + 3);
        constexpr size_t midstop_threshold = 11;
#endif
        if (m < midstop_threshold) {
            for (;;) {
                mask = 0;
                for (size_t i = start; i < end; i++) {
                    const FullVecUq first = bucket[i] + xVec;
                    const FullVecUi remapped = remap(first, first + FULL_VEC_64_COUNT, m);
                    mask |= powerOfTwo(remapped);
                }
#ifdef MORESTATS
                num_bij_evals[m] += m;
#endif
                if (horizontal_or(mask == found)) break;
                x += FULL_VEC_32_COUNT;
                xVec += FULL_VEC_32_COUNT;
            }
            const auto found_idx = horizontal_find_first(mask == found);
            assert(found_idx != -1);
            x += found_idx;
        } else {
            const uint32_t midstop = m / 2;
            uint32_t backlog_masks[FULL_VEC_32_COUNT];
            uint64_t backlog_x[FULL_VEC_32_COUNT];
            int backlog_size = 0;
            for (;;) {
                mask = 0;
                size_t i;
                for (i = start; i < start + midstop; i++) {
                    const FullVecUq first = bucket[i] + xVec;
                    const FullVecUi remapped = remap(first, first + FULL_VEC_64_COUNT, m);
                    mask |= powerOfTwo(remapped);
                }
#ifdef MORESTATS
                num_bij_evals[m] += midstop;
#endif
                bool brk = false;
#ifdef SIMDRS_512_BIT_POPCNT
                FullVecUi popcounts = _mm512_popcnt_epi32(mask);
                uint32_t bits = to_bits(popcounts == midstop);
                while (bits != 0) {
                    uint32_t j = __builtin_ctz(bits); // We don't care for MSVC here since it doesn't define __AVX512VPOPCNTDQ__
                    bits = clear_rho(bits);
                    uint32_t interm_mask = mask[j];
#else
                uint32_t intermediate[FULL_VEC_32_COUNT];
                mask.store(intermediate);
                for (uint32_t j = 0; j < FULL_VEC_32_COUNT; ++j) {
                    uint32_t interm_mask = intermediate[j];
                    if (nu(interm_mask) == (int)midstop) {
#endif
                        backlog_masks[backlog_size] = interm_mask;
                        backlog_x[backlog_size] = x + j;
                        if (++backlog_size == FULL_VEC_32_COUNT) {
                            FullVecUq first, second;
                            first.load(backlog_x);
                            second.load(backlog_x + FULL_VEC_64_COUNT);
                            FullVecUi backlog_mask;
                            backlog_mask.load(backlog_masks);
                            for (; i < end; i++) {
                                const FullVecUi remapped = remap(first + bucket[i], second + bucket[i], m);
                                backlog_mask |= powerOfTwo(remapped);
                            }
                            if (horizontal_or(backlog_mask == found)) {
                                brk = true;
                                const auto found_idx = horizontal_find_first(backlog_mask == found);
                                assert(found_idx != -1);
                                x = backlog_x[found_idx];
                                break;
                            }
                            backlog_size = 0;
                        }
#ifndef SIMDRS_512_BIT_POPCNT
                    }
#endif
                }
                if (brk) break;
                x += FULL_VEC_32_COUNT;
                xVec += FULL_VEC_32_COUNT;
            }
        }
#ifdef MORESTATS
        time_bij += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
        x -= SEED;
        const auto log2golomb = golomb_param(m);
        builder.appendFixed(x, log2golomb);
        unary.push_back(x >> log2golomb);
#ifdef MORESTATS
        bij_count[m]++;
        num_bij_trials[m] += x + 1;
        bij_unary += 1 + (x >> log2golomb);
        bij_fixed += log2golomb;

        min_bij_code = min(min_bij_code, x);
        max_bij_code = max(max_bij_code, x);
        sum_bij_codes += x;

        auto b = bij_memo_golomb[m];
        auto log2b = lambda(b);
        bij_unary_golomb += x / b + 1;
        bij_fixed_golomb += x % b < ((1 << (log2b + 1)) - b) ? log2b : log2b + 1;
#endif
    }

    // Computes and stores the splittings and bijections of a bucket.
    void recSplit(vector<uint64_t> &bucket, typename RiceBitVector<AT>::Builder &builder, vector<uint32_t> &unary) {
        const auto m = bucket.size();
        vector<uint64_t> temp(m);
        recSplit(bucket, temp, 0, bucket.size(), builder, unary, 0);
    }

    void recSplit(vector<uint64_t> &bucket, vector<uint64_t> &temp, size_t start, size_t m,
            typename RiceBitVector<AT>::Builder &builder, vector<uint32_t> &unary, const int level = 0) {
        assert(m > 1);
        if (m <= _leaf) {
            leafLevel(bucket, start, m, builder, unary, level);
        } else {
            [[maybe_unused]] uint64_t x;
            if (m > upper_aggr) { // fanout = 2
                const size_t split = ((uint16_t(m / 2 + upper_aggr - 1) / upper_aggr)) * upper_aggr;
                x = higherLevel(bucket, temp, start, m, split, builder, unary, level);
                recSplit(bucket, temp, start, split, builder, unary, level + 1);
                if (m - split > 1) recSplit(bucket, temp, start + split, m - split, builder, unary, level + 1);
#ifdef MORESTATS
                else
                    sum_depths += level;
#endif
            } else if (m > lower_aggr) { // 2nd aggregation level
                x = aggrLevel<upper_aggr / lower_aggr, lower_aggr, start_seed[NUM_START_SEEDS - 3]>(bucket, temp, start, m, builder, unary, level);
                size_t i;
                for (i = 0; i < m - lower_aggr; i += lower_aggr) {
                    recSplit(bucket, temp, start + i, lower_aggr, builder, unary, level + 1);
                }
                if (m - i > 1) recSplit(bucket, temp, start + i, m - i, builder, unary, level + 1);
#ifdef MORESTATS
                else
                    sum_depths += level;
#endif
            } else { // First aggregation level, m <= lower_aggr
                x = aggrLevel<lower_aggr / _leaf, _leaf, start_seed[NUM_START_SEEDS - 2]>(bucket, temp, start, m, builder, unary, level);
                size_t i;
                for (i = 0; i < m - _leaf; i += _leaf) {
                    leafLevel(bucket, start + i, _leaf, builder, unary, level + 1);
                }
                if (m - i > 1) leafLevel(bucket, start + i, m - i, builder, unary, level + 1);
#ifdef MORESTATS
                else
                    sum_depths += level;
#endif
            }

#ifdef MORESTATS
            ++split_count;
            num_split_trials += x + 1;
            double e_trials = 1;
            size_t aux = m;
            SplittingStrategy<LEAF_SIZE> strat(m);
            auto v = strat.begin();
            for (int i = 0; i < strat.fanout(); ++i, ++v) {
                e_trials *= pow((double)m / *v, *v);
                for (size_t j = *v; j > 0; --j, --aux) {
                    e_trials *= (double)j / aux;
                }
            }
            expected_split_trials += (size_t)e_trials;
            expected_split_evals += (size_t)e_trials * m;
            const auto log2golomb = golomb_param(m);
            split_unary += 1 + (x >> log2golomb);
            split_fixed += log2golomb;

            min_split_code = min(min_split_code, x);
            max_split_code = max(max_split_code, x);
            sum_split_codes += x;

            auto b = split_golomb_b<LEAF_SIZE>(m);
            auto log2b = lambda(b);
            split_unary_golomb += x / b + 1;
            split_fixed_golomb += x % b < ((1ULL << (log2b + 1)) - b) ? log2b : log2b + 1;
#endif
        }
    }

    void hash_gen(hash128_t *hashes) {
#ifdef MORESTATS
        time_bij = 0;
        memset(time_split, 0, sizeof time_split);
        split_unary = split_fixed = 0;
        bij_unary = bij_fixed = 0;
        min_split_code = 1ULL << 63;
        max_split_code = sum_split_codes = 0;
        min_bij_code = 1ULL << 63;
        max_bij_code = sum_bij_codes = 0;
        sum_depths = 0;
        minsize = this->keys_count;
        maxsize = 0;
        ub_split_bits = 0;
        ub_bij_bits = 0;
        ub_split_evals = 0;

        auto total_start_time = high_resolution_clock::now();
#endif

#ifndef __SIZEOF_INT128__
        if (this->keys_count > (1ULL << 32)) {
            fprintf(stderr, "For more than 2^32 keys, you need 128-bit integer support.\n");
            abort();
        }
#endif

        nbuckets = max(1, (keys_count + bucket_size - 1) / bucket_size);
        auto bucket_size_acc = vector<int64_t>(nbuckets + 1);
        auto bucket_pos_acc = vector<int64_t>(nbuckets + 1);
        ribbonInput.reserve(keys_count);

        sort_hash128_t(hashes, keys_count);
        typename RiceBitVector<AT>::Builder builder;

        bucket_size_acc[0] = bucket_pos_acc[0] = 0;
        for (size_t i = 0, last = 0; i < nbuckets; i++) {
            vector<uint64_t> bucket;
            for (; last < keys_count && hash128_to_bucket(hashes[last]) == i; last++) bucket.push_back(hashes[last].second);

            const size_t s = bucket.size();
            bucket_size_acc[i + 1] = bucket_size_acc[i] + s;
            if (bucket.size() > 1) {
                vector<uint32_t> unary;
                recSplit(bucket, builder, unary);
                builder.appendUnaryAll(unary);
            }
            bucket_pos_acc[i + 1] = builder.getBits();
        }
        builder.appendFixed(1, 1); // Sentinel (avoids checking for parts of size 1)
        descriptors = builder.build();
        ef = DoubleEF<AT>(vector<uint64_t>(bucket_size_acc.begin(), bucket_size_acc.end()), vector<uint64_t>(bucket_pos_acc.begin(), bucket_pos_acc.end()));

        // Begin: difference to RecSplit.
        ribbon = new SimpleRibbon<1>(ribbonInput);
        ribbonInput.clear();
        // End: difference to RecSplit.

#ifdef STATS
        // Evaluation purposes only
            double ef_sizes = (double)ef.bitCountCumKeys() / keys_count;
            double ef_bits = (double)ef.bitCountPosition() / keys_count;
            double rice_desc = (double)builder.getBits() / keys_count;
            double retrieval = 8.0 * (double)ribbon->size() / keys_count;
            printf("Elias-Fano cumul sizes:  %f bits/bucket\n", (double)ef.bitCountCumKeys() / nbuckets);
            printf("Elias-Fano cumul bits:   %f bits/bucket\n", (double)ef.bitCountPosition() / nbuckets);
            printf("Elias-Fano cumul sizes:  %f bits/key\n", ef_sizes);
            printf("Elias-Fano cumul bits:   %f bits/key\n", ef_bits);
            printf("Rice-Golomb descriptors: %f bits/key\n", rice_desc);
            printf("Retrieval:               %f bits/key\n", retrieval);
            printf("Total bits:              %f bits/key\n", ef_sizes + ef_bits + rice_desc + retrieval);

            printf("Total split bits        %16.3f\n", (double)split_fixed + split_unary);
            printf("Total bij bits:         %16.3f\n", (double)bij_fixed + bij_unary);

            printf("\n");
            printf("Bijections: %13.3f ms\n", time_bij * 1E-6);
            for (int i = 0; i < MAX_LEVEL_TIME; i++) {
                if (time_split[i] > 0) {
                    printf("Split level %d: %10.3f ms\n", i, time_split[i] * 1E-6);
                }
            }
#endif
    }

    static uint32_t remixAndRemap(uint64_t x, uint32_t n) {
        constexpr int masklen = 16;
        constexpr uint32_t mask = (uint32_t(1) << masklen) - 1;
        return ((remix32(uint32_t(x >> 32) ^ uint32_t(x)) & mask) * n) >> masklen;
    }

public:
    // TODO: why isn't this const?
    /** Returns the value associated with the given 128-bit hash.
     *
     * Note that this method is mainly useful for benchmarking.
     * @param hash a 128-bit hash.
     * @return the associated value.
     */
    size_t operator()(const hash128_t &hash) {
        const size_t bucket = hash128_to_bucket(hash);
        uint64_t cum_keys, cum_keys_next, bit_pos;
        ef.get(bucket, cum_keys, cum_keys_next, bit_pos);

        // Number of keys in this bucket
        size_t m = cum_keys_next - cum_keys;
        auto reader = descriptors.reader();
        reader.readReset(bit_pos, skip_bits(m));
        int level = 0;

        while (m > upper_aggr) { // fanout = 2
            const auto d = reader.readNext(golomb_param(m));
            const size_t hmod = remixAndRemap(hash.second + d + start_seed[level], m);

            const uint32_t split = ((uint16_t(m / 2 + upper_aggr - 1) / upper_aggr)) * upper_aggr;
            if (hmod < split) {
                m = split;
            } else {
                reader.skipSubtree(skip_nodes(split), skip_bits(split));
                m -= split;
                cum_keys += split;
            }
            level++;
        }
        if (m > lower_aggr) {
            const auto d = reader.readNext(golomb_param(m));
            const size_t hmod = remixAndRemap(hash.second + d + start_seed[NUM_START_SEEDS - 3], m);

            const int part = uint16_t(hmod) / lower_aggr;
            m = min(lower_aggr, m - part * lower_aggr);
            cum_keys += lower_aggr * part;
            if (part) reader.skipSubtree(skip_nodes(lower_aggr) * part, skip_bits(lower_aggr) * part);
        }

        if (m > _leaf) {
            const auto d = reader.readNext(golomb_param(m));
            const size_t hmod = remixAndRemap(hash.second + d + start_seed[NUM_START_SEEDS - 2], m);

            const int part = uint16_t(hmod) / _leaf;
            m = min(_leaf, m - part * _leaf);
            cum_keys += _leaf * part;
            if (part) reader.skipSubtree(part, skip_bits(_leaf) * part);
        }

        const auto b = reader.readNext(golomb_param(m));
        const uint64_t x = hash.second + b + start_seed[NUM_START_SEEDS - 1];
        return cum_keys + remixAndRemap(x, m);
    }

    /** Returns the value associated with the given key.
     *
     * @param key a key.
     * @return the associated value.
     */
    size_t operator()(const string &key) { return operator()(first_hash(key.c_str(), key.size())); }

    /** Returns an estimate of the size in bits of this structure. */
    size_t getBits() {
        return ef.bitCountCumKeys() + ef.bitCountPosition()
               + descriptors.getBits() + 8 * ribbon->size() + 8 * sizeof(SIMDShockHash);
    }

};

} // namespace shockhash