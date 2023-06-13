#pragma once

#include <vectorclass.h>

namespace shockhash {

using Vec4x64ui = Vec4uq;

/** David Stafford's (http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html)
* 13th variant of the 64-bit finalizer function in Austin Appleby's
* MurmurHash3 (https://github.com/aappleby/smhasher).
*/
Vec4x64ui inline remix(Vec4x64ui z) {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

Vec4x64ui remap32(Vec4x64ui remixed, uint32_t n) {
    constexpr uint32_t mask = (uint64_t(1) << 32) - 1;
    return ((remixed & mask) * n) >> 32;
}

Vec4x64ui powerOf2(Vec4x64ui x) {
    return _mm256_sllv_epi64(Vec4x64ui(1), x);
}

} // namespace shockhash
