#pragma once
/*
 * Based on RecSplit, Copyright (C) 2019-2020 Emmanuel Esposito and Sebastiano Vigna
 * Enhanced to use overloaded cuckoo hash tables in the leaves.
 * For tiny space usages (~1.6 bit/object), ShockHash is faster than RecSplit.
 */

#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <thread>
#include <condition_variable>
#include <sux/util/Vector.hpp>
#include <sux/function/DoubleEF.hpp>
#include <sux/function/RiceBitVector.hpp>
#include <sux/function/RecSplit.hpp>
#include <SimpleRibbon.h>
#include <Sorter.hpp>
#include "ShockHash.h"
#include "ShockHash2-internal.h"
#include "RiceBitVector.h"

namespace shockhash {
static const int MAX_LEAF_SIZE2 = 138;

// Optimal Golomb-Rice parameters for leaves.
static constexpr uint8_t bij_memo2[MAX_LEAF_SIZE2 + 1] = {
        0,  0,  0,  0,  0,  0,  1,  1,  2,  2, // 0..9
        2,  3,  3,  4,  4,  4,  5,  5,  6,  6, // 10..19
        7,  7,  8,  8,  9,  9,  9, 10, 10, 10, // 20..29
        12, 12, 13, 13, 14, 14, 14, 15, 15, 16, // 30..39
        16, 17, 17, 17, 18, 18, 19, 19, 20, 20, // 40..49
        20, 21, 21, 21, 22, 22, 23, 23, 24, 24, // 50..59
        24, 24, 25, 25, 25, 26, 26, 26, 27, 27, // 60..69
        28, 28, 28, 29, 30, 30, 30, 31, 32, 31, // 70..79
        32, 32, 33, 33, 33, 34, 35, 35, 36, 36, // 80..89
        36, 36, 37, 37, 38, 38, 39, 39, 40, 40, // 90..99
        41, 41, 43, 43, 44, 44, 44, 44, 45, 45, // 100..109
        46, 46, 46, 47, 47, 47, 48, 48, 49, 49, // 110..119
        50, 50, 50, 51, 51, 52, 52, 53, 53, 54, // 120..129
        54, 55, 55, 56, 56, 57, 58, 58, 58, // 130..138
};

template <size_t LEAF_SIZE> class SplittingStrategy2 {
    public:
        static constexpr size_t _leaf = LEAF_SIZE;
        static_assert(_leaf >= 1);
        static_assert(_leaf <= MAX_LEAF_SIZE2);
        static constexpr size_t lower_aggr = _leaf * 4;
        static constexpr size_t upper_aggr = lower_aggr * 3;
};

// Generates the precomputed table of 32-bit values holding the Golomb-Rice code
// of a splitting (upper 5 bits), the number of nodes in the associated subtree
// (following 11 bits) and the sum of the Golomb-Rice codelengths in the same
// subtree (lower 16 bits).

template <size_t LEAF_SIZE> static constexpr void _fill_golomb_rice2(const size_t m, array<uint64_t, MAX_BUCKET_SIZE> *memo) {
    array<long, MAX_FANOUT> k{0};

    constexpr size_t lower_aggr = SplittingStrategy2<LEAF_SIZE>::lower_aggr;
    constexpr size_t upper_aggr = SplittingStrategy2<LEAF_SIZE>::upper_aggr;

    size_t fanout = 0, unit = 0;
    if (m > upper_aggr) { // High-level aggregation (fanout 2)
        unit = upper_aggr * (uint16_t(m / 2 + upper_aggr - 1) / upper_aggr);
        fanout = 2;
    } else if (m > lower_aggr) { // Second-level aggregation
        unit = lower_aggr;
        fanout = uint16_t(m + lower_aggr - 1) / lower_aggr;
    } else { // First-level aggregation
        unit = LEAF_SIZE;
        fanout = uint16_t(m + LEAF_SIZE - 1) / LEAF_SIZE;
    }

    k[fanout - 1] = m;
    for (size_t i = 0; i < fanout - 1; ++i) {
        k[i] = unit;
        k[fanout - 1] -= k[i];
    }

    double sqrt_prod = 1;
    for (size_t i = 0; i < fanout; ++i) sqrt_prod *= sqrt(k[i]);

    const double p = sqrt(m) / (pow(2 * M_PI, (fanout - 1.) / 2) * sqrt_prod);
    uint64_t golomb_rice_length = ceil(log2(-log((sqrt(5) + 1) / 2) / log1p(-p))); // log2 Golomb modulus

    assert(golomb_rice_length <= 0x1F); // Golomb-Rice code, stored in the 5 upper bits
    assert((golomb_rice_length << 27) >> 27 == golomb_rice_length);
    (*memo)[m] = golomb_rice_length << 27;
    for (size_t i = 0; i < fanout; ++i) golomb_rice_length += (*memo)[k[i]] & 0xFFFF;
    assert(golomb_rice_length <= 0xFFFF); // Sum of Golomb-Rice codeslengths in the subtree, stored in the lower 16 bits
    (*memo)[m] |= golomb_rice_length;

    uint32_t nodes = 1;
    for (size_t i = 0; i < fanout; ++i) nodes += ((*memo)[k[i]] >> 16) & 0x7FF;
    assert(LEAF_SIZE < 3 || nodes <= 0x7FF); // Number of nodes in the subtree, stored in the middle 11 bits
    (*memo)[m] |= nodes << 16;
}

template <size_t LEAF_SIZE> static constexpr array<uint64_t, MAX_BUCKET_SIZE> fill_golomb_rice2() {
    array<uint64_t, MAX_BUCKET_SIZE> memo{0};
    size_t s = 0;
    for (; s <= LEAF_SIZE; ++s) {
        memo[s] = uint64_t(bij_memo2[s]) << 27 | (s > 1) << 16 | bij_memo2[s];
        assert(memo[s] >> 27 == bij_memo2[s]);
    }
    for (; s < MAX_BUCKET_SIZE; ++s) _fill_golomb_rice2<LEAF_SIZE>(s, &memo);
    return memo;
}

template <size_t LEAF_SIZE>
class ShockHash2 {
        static_assert(LEAF_SIZE <= MAX_LEAF_SIZE2);
        static constexpr AllocType AT = sux::util::AllocType::MALLOC;
        static constexpr size_t _leaf = LEAF_SIZE;
        static constexpr size_t lower_aggr = SplittingStrategy2<LEAF_SIZE>::lower_aggr;
        static constexpr size_t upper_aggr = SplittingStrategy2<LEAF_SIZE>::upper_aggr;

        // For each bucket size, the Golomb-Rice parameter (upper 8 bits) and the number of bits to
        // skip in the fixed part of the tree (lower 24 bits).
        static constexpr array<uint64_t, MAX_BUCKET_SIZE> memo = fill_golomb_rice2<LEAF_SIZE>();

        size_t bucket_size;
        size_t nbuckets;
        size_t keys_count;
        RiceBitVector<AT> descriptors;
        DoubleEF<AT> ef;
        using Ribbon = SimpleRibbon<1, (_leaf > 24) ? 128 : 64>;
        Ribbon *ribbon = nullptr;
        std::vector<std::pair<uint64_t, uint8_t>> ribbonInput;
        std::mutex ribbonInputMtx;

    public:
        ShockHash2() {}


        ShockHash2(const vector<string> &keys, const size_t bucket_size, size_t num_threads = 1) {
            this->bucket_size = bucket_size;
            this->keys_count = keys.size();
            hash128_t *h = (hash128_t *)malloc(this->keys_count * sizeof(hash128_t));
            if (num_threads == 1) {
                for (size_t i = 0; i < this->keys_count; ++i) {
                    h[i] = first_hash(keys[i].c_str(), keys[i].size());
                }
            } else {
                size_t keysPerThread = this->keys_count / num_threads + 1;
                std::vector<std::thread> threads;
                for (size_t thread = 0; thread < num_threads; thread++) {
                    threads.emplace_back([&, thread] {
                        size_t from = thread * keysPerThread;
                        size_t to = std::min(this->keys_count, (thread + 1) * keysPerThread);
                        for (size_t i = from; i < to; ++i) {
                            h[i] = first_hash(keys[i].c_str(), keys[i].size());
                        }
                    });
                }
                for (std::thread &t : threads) {
                    t.join();
                }
            }
            hash_gen(h, num_threads);
            free(h);
        }

        ShockHash2(vector<hash128_t> &keys, const size_t bucket_size, size_t num_threads = 1) {
            this->bucket_size = bucket_size;
            this->keys_count = keys.size();
            hash_gen(&keys[0], num_threads);
        }

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
                const size_t hmod = sux::remap16(sux::function::remix(hash.second + d + start_seed[level]), m);

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
                const size_t hmod = sux::remap16(sux::function::remix(hash.second + d + start_seed[level]), m);

                const int part = uint16_t(hmod) / lower_aggr;
                m = min(lower_aggr, m - part * lower_aggr);
                cum_keys += lower_aggr * part;
                if (part) reader.skipSubtree(skip_nodes(lower_aggr) * part, skip_bits(lower_aggr) * part);
                level++;
            }

            if (m > _leaf) {
                const auto d = reader.readNext(golomb_param(m));
                const size_t hmod = sux::remap16(sux::function::remix(hash.second + d + start_seed[level]), m);

                const int part = uint16_t(hmod) / _leaf;
                m = min(_leaf, m - part * _leaf);
                cum_keys += _leaf * part;
                if (part) reader.skipSubtree(part, skip_bits(_leaf) * part);
                level++;
            }

            const auto b = reader.readNext(golomb_param(m));

            // Begin: difference to RecSplit.
            return cum_keys + dispatchLeafSizeQ<LEAF_SIZE>(m, b, hash.second, ribbon->retrieve(hash.second));
            // End: difference to RecSplit.
        }

        template <size_t I>
        inline size_t dispatchLeafSizeQ(size_t param, size_t seed, uint64_t key, size_t retrieved) {
            if constexpr (I <= 1) {
                return 0;
            } else if (I == param) {
                using SH = std::conditional_t<(I >= 10),
                        BijectionsShockHash2<I, true, QuadSplitCandidateFinderTree>,
                        BijectionsShockHash2<I, true, BasicSeedCandidateFinder>>;
                return SH::hash(seed, key, retrieved);
            } else {
                return dispatchLeafSizeQ<I - 1>(param, seed, key, retrieved);
            }
        }

        /** Returns the value associated with the given key.
         *
         * @param key a key.
         * @return the associated value.
         */
        size_t operator()(const string &key) { return operator()(first_hash(key.c_str(), key.size())); }

        /** Returns the number of keys used to build this RecSplit instance. */
        inline size_t size() { return this->keys_count; }

        /** Returns an estimate of the size in bits of this structure. */
        size_t getBits() {
            return ef.bitCountCumKeys() + ef.bitCountPosition()
                    + descriptors.getBits() + 8 * ribbon->size() + 8 * sizeof(ShockHash2);
        }

        void printBits() {
            std::cout<<"EF 1:   "<<(double)ef.bitCountCumKeys()/keys_count<<std::endl;
            std::cout<<"EF 2:   "<<(double)ef.bitCountPosition()/keys_count<<std::endl;
            std::cout<<"trees:  "<<(double)descriptors.getBits()/keys_count<<std::endl;
            std::cout<<"ribbon: "<<(double)(8 * ribbon->size())/keys_count<<std::endl;
        }

    private:
        // Maps a 128-bit to a bucket using the first 64-bit half.
        inline uint64_t hash128_to_bucket(const hash128_t &hash) const { return remap128(hash.first, nbuckets); }

        template <size_t I>
        inline size_t dispatchLeafSize(size_t param, std::vector<uint64_t> &leafKeys) {
            if constexpr (I <= 1) {
                return 0;
            } else if (I == param) {
                using SH = std::conditional_t<(I >= 10),
                        BijectionsShockHash2<I, true, QuadSplitCandidateFinderTree>,
                        BijectionsShockHash2<I, true, BasicSeedCandidateFinder>>;
                size_t x = SH::findSeed(leafKeys);
                {
                    std::lock_guard<std::mutex> guard(ribbonInputMtx);
                    SH::constructRetrieval(leafKeys, x, ribbonInput);
                }
                #ifndef NDEBUG
                    SH::verify(x, leafKeys, ribbonInput);
                #endif
                return x;
            } else {
                return dispatchLeafSize<I - 1>(param, leafKeys);
            }
        }

        void recSplit(vector<uint64_t> &bucket, vector<uint64_t> &temp, size_t start, size_t end,
                      typename RiceBitVector<AT>::Builder &builder, vector<uint32_t> &unary, const int level,
                      TinyBinaryCuckooHashTable &tinyBinaryCuckooHashTable) {
            const auto m = end - start;
            assert(m > 1);
            uint64_t x = start_seed[level];

            if (m <= _leaf) {
#ifdef STATS
                auto start_time = high_resolution_clock::now();
#endif
                // Begin: difference to RecSplit.
                std::vector<uint64_t> leafKeys(bucket.begin() + start, bucket.begin() + end);
                x = dispatchLeafSize<LEAF_SIZE>(m, leafKeys);
                // End: difference to RecSplit.

                const auto log2golomb = golomb_param(m);
                assert(log2golomb > 0);
                builder.appendFixed(x, log2golomb);
                unary.push_back(x >> log2golomb);

#ifdef STATS
                bij_unary += 1 + (x >> log2golomb);
                bij_fixed += log2golomb;
                time_bij += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
                bij_opt += log2(x + 1);
#endif
            } else {
#ifdef STATS
                auto start_time = high_resolution_clock::now();
#endif
                if (m > upper_aggr) { // fanout = 2
                    const size_t split = ((uint16_t(m / 2 + upper_aggr - 1) / upper_aggr)) * upper_aggr;

                    size_t count[2];
                    for (;;) {
                        count[0] = 0;
                        for (size_t i = start; i < end; i++) {
                            count[remap16(sux::function::remix(bucket[i] + x), m) >= split]++;
                        }
                        if (count[0] == split) break;
                        x++;
                    }

                    count[0] = 0;
                    count[1] = split;
                    for (size_t i = start; i < end; i++) {
                        temp[count[remap16(sux::function::remix(bucket[i] + x), m) >= split]++] = bucket[i];
                    }
                    copy(&temp[0], &temp[m], &bucket[start]);
                    x -= start_seed[level];
                    const auto log2golomb = golomb_param(m);
                    builder.appendFixed(x, log2golomb);
                    unary.push_back(x >> log2golomb);

#ifdef STATS
                    time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
                    split_opt += log2(x + 1);
#endif
                    recSplit(bucket, temp, start, start + split, builder, unary, level + 1, tinyBinaryCuckooHashTable);
                    if (m - split > 1) recSplit(bucket, temp, start + split, end, builder, unary, level + 1, tinyBinaryCuckooHashTable);
                } else if (m > lower_aggr) { // 2nd aggregation level
                    const size_t fanout = uint16_t(m + lower_aggr - 1) / lower_aggr;
                    size_t count[fanout]; // Note that we never read count[fanout-1]
                    for (;;) {
                        memset(count, 0, sizeof count - sizeof *count);
                        for (size_t i = start; i < end; i++) {
                            count[uint16_t(remap16(sux::function::remix(bucket[i] + x), m)) / lower_aggr]++;
                        }
                        size_t broken = 0;
                        for (size_t i = 0; i < fanout - 1; i++) broken |= count[i] - lower_aggr;
                        if (!broken) break;
                        x++;
                    }

                    for (size_t i = 0, c = 0; i < fanout; i++, c += lower_aggr) count[i] = c;
                    for (size_t i = start; i < end; i++) {
                        temp[count[uint16_t(sux::function::remap16(sux::function::remix(bucket[i] + x), m)) / lower_aggr]++] = bucket[i];
                    }
                    copy(&temp[0], &temp[m], &bucket[start]);

                    x -= start_seed[level];
                    const auto log2golomb = golomb_param(m);
                    builder.appendFixed(x, log2golomb);
                    unary.push_back(x >> log2golomb);

#ifdef STATS
                    time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
                    split_opt += log2(x + 1);
#endif
                    size_t i;
                    for (i = 0; i < m - lower_aggr; i += lower_aggr) {
                        recSplit(bucket, temp, start + i, start + i + lower_aggr, builder, unary, level + 1, tinyBinaryCuckooHashTable);
                    }
                    if (m - i > 1) recSplit(bucket, temp, start + i, end, builder, unary, level + 1, tinyBinaryCuckooHashTable);
                } else { // First aggregation level, m <= lower_aggr
                    const size_t fanout = uint16_t(m + _leaf - 1) / _leaf;
                    size_t count[fanout]; // Note that we never read count[fanout-1]
                    for (;;) {
                        memset(count, 0, sizeof count - sizeof *count);
                        for (size_t i = start; i < end; i++) {
                            count[uint16_t(remap16(sux::function::remix(bucket[i] + x), m)) / _leaf]++;
                        }
                        size_t broken = 0;
                        for (size_t i = 0; i < fanout - 1; i++) broken |= count[i] - _leaf;
                        if (!broken) break;
                        x++;
                    }
                    for (size_t i = 0, c = 0; i < fanout; i++, c += _leaf) count[i] = c;
                    for (size_t i = start; i < end; i++) {
                        temp[count[uint16_t(remap16(sux::function::remix(bucket[i] + x), m)) / _leaf]++] = bucket[i];
                    }
                    copy(&temp[0], &temp[m], &bucket[start]);

                    x -= start_seed[level];
                    const auto log2golomb = golomb_param(m);
                    builder.appendFixed(x, log2golomb);
                    unary.push_back(x >> log2golomb);

#ifdef STATS
                    time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
                    split_opt += log2(x + 1);
#endif
                    size_t i;
                    for (i = 0; i < m - _leaf; i += _leaf) {
                        recSplit(bucket, temp, start + i, start + i + _leaf, builder, unary, level + 1, tinyBinaryCuckooHashTable);
                    }
                    if (m - i > 1) recSplit(bucket, temp, start + i, end, builder, unary, level + 1, tinyBinaryCuckooHashTable);
                }
#ifdef STATS
                const auto log2golomb = golomb_param(m);
                split_unary += 1 + (x >> log2golomb);
                split_fixed += log2golomb;
#endif
            }
        }

        void compute_thread(int tid, int num_threads, mutex &mtx, std::condition_variable &condition,
                            vector<uint64_t> &bucket_size_acc, vector<uint64_t> &bucket_pos_acc,
                            vector<uint64_t> &sorted_keys, int &next_thread_to_append_builder,
                            typename shockhash::RiceBitVector<AT>::Builder &builder) {
            typename RiceBitVector<AT>::Builder local_builder;
            TinyBinaryCuckooHashTable tinyBinaryCuckooHashTable(LEAF_SIZE);
            vector<uint32_t> unary;
            vector<uint64_t> temp(MAX_BUCKET_SIZE);
            size_t begin = tid * this->nbuckets / num_threads;
            size_t end = std::min(this->nbuckets, (tid + 1) * this->nbuckets / num_threads);
            if (tid == num_threads - 1) {
                end = this->nbuckets;
            }
            for (size_t i = begin; i < end; ++i) {
                const size_t s = bucket_size_acc[i + 1] - bucket_size_acc[i];
                if (s > 1) {
                    recSplit(sorted_keys, temp, bucket_size_acc[i], bucket_size_acc[i + 1], local_builder,
                             unary, 0, tinyBinaryCuckooHashTable);
                    local_builder.appendUnaryAll(unary);
                    unary.clear();
                }
                bucket_pos_acc[i + 1] = local_builder.getBits();
            }
            if (tid == 0) {
                builder = std::move(local_builder);
                lock_guard<mutex> lock(mtx);
                next_thread_to_append_builder = 1;
                condition.notify_all();
            } else {
                uint64_t prev_bucket_pos;
                {
                    unique_lock<mutex> lock(mtx);
                    condition.wait(lock, [&] { return next_thread_to_append_builder == tid; });
                    prev_bucket_pos = builder.getBits();
                    builder.appendRiceBitVector(local_builder);
                    next_thread_to_append_builder = tid + 1;
                    condition.notify_all();
                }
                for (size_t i = begin + 1; i < end + 1; ++i) {
                    bucket_pos_acc[i] += prev_bucket_pos;
                }
            }
        }

        void hash_gen(hash128_t *hashes, int num_threads) {
#ifdef STATS
            split_unary = split_fixed = 0;
            bij_unary = bij_fixed = 0;
            time_bij = 0;
            memset(time_split, 0, sizeof time_split);
#endif

#ifndef __SIZEOF_INT128__
            if (keys_count > (1ULL << 32)) {
			fprintf(stderr, "For more than 2^32 keys, you need 128-bit integer support.\n");
			abort();
		}
#endif
            nbuckets = max(1, (keys_count + bucket_size - 1) / bucket_size);
            auto bucket_size_acc = std::vector<uint64_t>(nbuckets + 1);
            auto bucket_pos_acc = std::vector<uint64_t>(nbuckets + 1);
            auto sorted_keys = vector<uint64_t>(keys_count);
            ribbonInput.reserve(keys_count);

            parallelPartition(hashes, sorted_keys, bucket_size_acc, num_threads, keys_count, nbuckets);
            typename RiceBitVector<AT>::Builder builder;

            vector<std::thread> threads;
            threads.reserve(num_threads);
            mutex mtx;
            std::condition_variable condition;
            int next_thread_to_append_builder = 0;
            bucket_pos_acc[0] = 0;
            if (num_threads == 1) {
                compute_thread(0, num_threads, mtx, condition,
                               bucket_size_acc, bucket_pos_acc, sorted_keys,
                               next_thread_to_append_builder, builder);
            } else {
                for (int tid = 0; tid < num_threads; ++tid) {
                    threads.emplace_back([&, tid] {
                        compute_thread(tid, num_threads, mtx, condition,
                                       bucket_size_acc, bucket_pos_acc, sorted_keys,
                                       next_thread_to_append_builder, builder);
                    });
                }
                for (auto &thread: threads) {
                    thread.join();
                }
            }
            builder.appendFixed(1, 1); // Sentinel (avoids checking for parts of size 1)
            descriptors = builder.build();
            ef = DoubleEF<AT>(vector<uint64_t>(bucket_size_acc.begin(), bucket_size_acc.end()), vector<uint64_t>(bucket_pos_acc.begin(), bucket_pos_acc.end()));

            // Begin: difference to RecSplit.
            ribbon = new Ribbon(ribbonInput);
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

            printf("Split bits:       %16.3f\n", ((double)split_fixed + split_unary) / keys_count);
            printf("Split bits opt:   %16.3f\n", split_opt / keys_count);
            printf("Bij bits:         %16.3f\n", ((double)bij_fixed + bij_unary) / keys_count);
            printf("Bij bits opt:     %16.3f\n", bij_opt / keys_count);

            printf("\n");
            printf("Bijections: %13.3f ms\n", time_bij * 1E-6);
            //for (size_t i = 0; i <= LEAF_SIZE; i++) {
            //    printf("Bijections of size %d:    %d\n", i, bij_count[i]);
            //}
            for (int i = 0; i < MAX_LEVEL_TIME; i++) {
                if (time_split[i] > 0) {
                    printf("Split level %d: %10.3f ms\n", i, time_split[i] * 1E-6);
                }
            }
#endif
        }
};

} // namespace shockhash
