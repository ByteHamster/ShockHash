#pragma once

#include <vectorclass.h>

namespace shockhash {

#if INSTRSET >= 9
    #define SHOCKHASH_SIMD_512_BIT
#elif INSTRSET == 8
    #define SHOCKHASH_SIMD_256_BIT
#else
    #pragma message("ShockHash was compiled without AVX512 and AVX2 support => suboptimal performance")
#endif

#ifdef SHOCKHASH_SIMD_512_BIT
    using FullVecUi = Vec16ui; // 16x32ui
    using FullVecUq = Vec8uq;  // 8x64ui
    using FullVecIb = Vec16ib;
    using FullVecQ = Vec8q;
    using FullVecC = Vec64c;
    static constexpr uint32_t FULL_VEC_32_COUNT = 16;
    static const FullVecUq FULL_VEC_COUNTING = FullVecUq(0, 1, 2, 3, 4, 5, 6, 7);
#else
    using FullVecUi = Vec8ui; // 8x32ui
    using FullVecUq = Vec4uq; // 4x64ui
    using FullVecIb = Vec8ib;
    using FullVecQ = Vec4q;
    using FullVecC = Vec32c;
    static constexpr uint32_t FULL_VEC_32_COUNT = 8;
    static const FullVecUq FULL_VEC_COUNTING = FullVecUq(0, 1, 2, 3);
#endif

static constexpr uint32_t FULL_VEC_64_COUNT = FULL_VEC_32_COUNT / 2;
static const FullVecUq FULL_VEC_ALL_ONE = FullVecUq(1);
static const FullVecUq FULL_VEC_ALL_ZERO = FullVecUq(0);
static const FullVecUq FULL_VEC_ALL_64_COUNT = FullVecUq(FULL_VEC_64_COUNT);

/** David Stafford's (http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html)
* 13th variant of the 64-bit finalizer function in Austin Appleby's
* MurmurHash3 (https://github.com/aappleby/smhasher).
*/
FullVecUq inline remixV(FullVecUq z) {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

FullVecUq remap32V(FullVecUq remixed, uint32_t n) {
    constexpr uint32_t mask = (uint64_t(1) << 32) - 1;
    return ((remixed & mask) * n) >> 32;
}

static FullVecUq shiftLeftV(FullVecUq x, FullVecUq y) {
    #ifdef SHOCKHASH_SIMD_512_BIT
        return _mm512_sllv_epi64(x, y);
    #else
        return _mm256_sllv_epi64(x, y);
    #endif
}

static FullVecUq shiftRightV(FullVecUq x, FullVecUq y) {
    #ifdef SHOCKHASH_SIMD_512_BIT
        return _mm512_srlv_epi64(x, y);
    #else
        return _mm256_srlv_epi64(x, y);
    #endif
}

static FullVecUq powerOfTwo(FullVecUq x) {
    return shiftLeftV(FULL_VEC_ALL_ONE, x);
}

} // namespace shockhash
