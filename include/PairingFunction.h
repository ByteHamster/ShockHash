#pragma once

#include <cmath>

namespace shockhash {
/**
 * Cantor Pairing Function: https://en.wikipedia.org/wiki/Pairing_function
 * Assigns numbers diagonally
 *   y
 * 3 | 9
 * 2 | 5 8
 * 1 | 2 4 7
 * 0 | 0 1 3 6
 *    --------- x
 *     0 1 2 3
 */
inline size_t pairCantor(size_t x, size_t y) {
    return (x + y) * (x + y + 1) / 2 + y;
}

inline std::pair<size_t, size_t> unpairCantor(size_t z) {
    size_t w = std::floor((sqrt(8.0 * (double)z + 1) - 1) / 2);
    size_t t = (w * w + w) / 2;
    size_t y = z - t;
    size_t x = w - y;
    return std::make_pair(x, y);
}

inline bool isFloatAccurateCantor(size_t x, size_t y) {
    auto [restoredX, restoredY] = unpairCantor(pairCantor(x, y));
    return x == restoredX && y == restoredY;
}

/**
 * Our own hand-crafted pairing function.
 * Assumes x > y and does not use diagonals.
 * Not using diagonals helps with lazy evaluation if x and y come from the same set.
 *
 *     0 1 2 3 y
 * 0 | -
 * 1 | 0
 * 2 | 1 2
 * 3 | 3 4 5
 * 4 | 6 7 8 9
 * x
 */
inline size_t pairTriangular(size_t x, size_t y) {
    assert(x > y);
    return (x - 1) * x / 2 + y;
}

inline std::pair<size_t, size_t> unpairTriangular(size_t z) {
    uint64_t x = floor(0.5 + sqrt(0.25 + 2.0 * (double)z));
    uint64_t y = z - ((x - 1) * x) / 2;
    return std::make_pair(x, y);
}

inline bool isFloatAccurateTriangular(size_t x, size_t y) {
    auto [restoredX, restoredY] = unpairTriangular(pairTriangular(x, y));
    return x == restoredX && y == restoredY;
}

/**
 * 2006, Szudzik: Elegant Pair
 * Assigns numbers along the edges of squares.
 * Helps with lazy evaluation if always calculating the same new value for both x and y.
 *   y
 * 3 |  9 10 11 15
 * 2 |  4  5  8 14
 * 1 |  1  3  7 13
 * 0 |  0  2  6 12
 *    ------------ x
 *      0  1  2  3
 */
inline size_t pairElegant(size_t x, size_t y) {
    if (x < y) {
        return y * y + x;
    } else {
        return x * x + x + y;
    }
}

inline std::pair<size_t, size_t> unpairElegant(size_t z) {
    size_t floorSqrt = floor(sqrt((double) z));
    size_t reSquared = floorSqrt * floorSqrt;
    if (z - reSquared < floorSqrt) {
        return std::make_pair(z - reSquared, floorSqrt);
    } else {
        return std::make_pair(floorSqrt, z - reSquared - floorSqrt);
    }
}

inline bool isFloatAccurateElegant(size_t x, size_t y) {
    auto [restoredX, restoredY] = unpairElegant(pairElegant(x, y));
    return x == restoredX && y == restoredY;
}
} // namespace shockhash
