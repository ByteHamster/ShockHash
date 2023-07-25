#pragma once

#include <vectorclass.h>

namespace shockhash {

#if INSTRSET >= 9
#define SIMDRS_512_BIT
constexpr bool NO_AVX = false;

#ifdef __AVX512VPOPCNTDQ__
#define SIMDRS_512_BIT_POPCNT
#pragma message("AVX512 with popcount")
#else
#pragma message("AVX512")
#endif

#elif INSTRSET == 8
constexpr bool NO_AVX = false;
#else
constexpr bool NO_AVX = true;
#pragma message("SIMDRecSplit was compiled without AVX512 and AVX2 support => suboptimal performance")
#endif

#ifdef SIMDRS_512_BIT
using FullVecUi = Vec16ui;
using FullVecUq = Vec8uq;
using FullVecIb = Vec16ib;
using FullVecQ = Vec8q;
using FullVecC = Vec64c;
constexpr uint32_t FULL_VEC_32_COUNT = 16;
#else
using FullVecUi = Vec8ui;
using FullVecUq = Vec4uq;
using FullVecIb = Vec8ib;
using FullVecQ = Vec4q;
using FullVecC = Vec32c;
constexpr uint32_t FULL_VEC_32_COUNT = 8;
#endif
constexpr uint32_t FULL_VEC_64_COUNT = FULL_VEC_32_COUNT / 2;

using Vec4x64ui = Vec4uq;

static const Vec4x64ui VEC_0123 = Vec4x64ui(0, 1, 2, 3);
static const Vec4x64ui VEC_1111 = Vec4x64ui(1);

/** David Stafford's (http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html)
* 13th variant of the 64-bit finalizer function in Austin Appleby's
* MurmurHash3 (https://github.com/aappleby/smhasher).
*/
Vec4x64ui inline remixV(Vec4x64ui z) {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

Vec4x64ui remap32V(Vec4x64ui remixed, uint32_t n) {
    constexpr uint32_t mask = (uint64_t(1) << 32) - 1;
    return ((remixed & mask) * n) >> 32;
}

static Vec4x64ui shiftLeftV(Vec4x64ui x, Vec4x64ui y) {
    return _mm256_sllv_epi64(x, y);
}

static Vec4x64ui shiftRightV(Vec4x64ui x, Vec4x64ui y) {
    return _mm256_srlv_epi64(x, y);
}

static Vec4x64ui powerOfTwo(Vec4x64ui x) {
    return shiftLeftV(VEC_1111, x);
}

} // namespace shockhash
