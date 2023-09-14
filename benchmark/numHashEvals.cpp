#include <iostream>
#include <vector>
#include <cmath>
#include "../include/TinyBinaryCuckooHashTable.h"
#include <XorShift64.h>
#include <chrono>
#include <set>
#include <unordered_set>

static constexpr uint64_t rotate(size_t l, uint64_t val, uint32_t x) {
    return ((val << x) | (val >> (l - x))) & ((1ul << l) - 1);
}

void testRotationFitting(size_t l) {
    size_t numIterations = std::max(2ul, (size_t) 7e7 / (1<<l));
    size_t totalTries = 0;
    size_t hfEvals = 0;
    for (size_t iteration = 0; iteration < numIterations; iteration++) {
        std::vector<shockhash::HashedKey> keys;
        for (size_t i = 0; i < l; i++) {
            keys.emplace_back(std::to_string(i) + " " + std::to_string(iteration));
        }
        while (true) {
            uint64_t taken1 = 0;
            uint64_t taken2 = 0;
            for (size_t i = 0; i < l; i++) {
                size_t pos = keys[i].hash(totalTries, l);
                hfEvals++;
                if (keys[i].mhc % 2 == 0) {
                    taken1 |= (1<<pos);
                } else {
                    taken2 |= (1<<pos);
                }
            }
            bool canBeRotated = false;
            for (size_t r = 0; r < l; r++) {
                if ((taken1 | rotate(l, taken2, r)) == (1<<l) - 1) {
                    canBeRotated = true;
                    break;
                }
            }
            if (canBeRotated) {
                break;
            }
            totalTries++;
        }
    }
    std::cout<<"RESULT"
             <<" method=rotations"
             <<" l="<<l
             <<" hfEvals="<<(double)hfEvals/(double)numIterations
             <<" tries="<<(double)totalTries / (double)numIterations
             <<" iterations="<<numIterations
             <<" spaceEstimate="<<log2((double)totalTries / (double)numIterations * l) / l
             <<std::endl;

}

void testBruteForce(size_t l) {
    size_t numIterations = std::max(2ul, (size_t) 4e7 / (1<<l));
    size_t totalTries = 0;
    size_t hfEvals = 0;
    for (size_t iteration = 0; iteration < numIterations; iteration++) {
        std::vector<shockhash::HashedKey> keys;
        for (size_t i = 0; i < l; i++) {
            keys.emplace_back(std::to_string(i) + " " + std::to_string(iteration));
        }
        while (true) {
            uint64_t taken = 0;
            size_t i = 0;
            for (; i < l; i++) {
                size_t pos = keys[i].hash(totalTries, l);
                hfEvals++;
                if (taken & (1<<pos)) {
                    break;
                }
                taken |= (1<<pos);
            }
            if (i == l) {
                break;
            }
            totalTries++;
        }
    }
    std::cout<<"RESULT"
             <<" method=bruteforce"
             <<" l="<<l
             <<" hfEvals="<<(double)hfEvals/(double)numIterations
             <<" tries="<<(double)totalTries / (double)numIterations
             <<" iterations="<<numIterations
             <<" spaceEstimate="<<log2((double)totalTries / (double)numIterations) / l
             <<std::endl;
}

void testShockHash(size_t l) {
    size_t numIterations = 20000;
    if (l > 50) {
        numIterations = 100;
    } else if (l > 42) {
        numIterations = 1000;
    } else if (l > 30) {
        numIterations = 4000;
    }
    size_t totalTries = 0;
    size_t hfEvals = 0;
    size_t totalPseudotrees = 0;
    size_t totalOrientations = 0;
    size_t totalWith1Component = 0;
    for (size_t iteration = 0; iteration < numIterations; iteration++) {
        shockhash::TinyBinaryCuckooHashTable table(l);
        for (size_t i = 0; i < l; i++) {
            table.prepare(shockhash::HashedKey(std::to_string(i) + " " + std::to_string(iteration)));
        }
        while (!table.testPassesFilter(totalTries) || !table.construct(totalTries)) { // totalTries is the (unique) seed here
            totalTries++;
            hfEvals += 2 * l;
        }
        // Find number of pseudotrees
        shockhash::UnionFind unionFind(l);
        for (size_t i = 0; i < l; i++) {
            shockhash::TinyBinaryCuckooHashTable::CandidateCells candidateCells
                    = shockhash::TinyBinaryCuckooHashTable::getCandidateCells(table.heap[i].hash, totalTries, l);
            bool res = unionFind.unionIsStillPseudoforrest(candidateCells.cell1, candidateCells.cell2);
            if (!res) {
                std::cerr << "Something went wrong" << std::endl;
                exit(1);
            }
        }
        std::unordered_set<uint64_t> representatives;
        for (size_t i = 0; i < l; i++) {
            size_t repr = unionFind.findRepresentative(i);
            representatives.insert(repr);
        }
        totalPseudotrees += representatives.size();
        totalOrientations += 1ul<<representatives.size();
        if (representatives.size() == 1) {
            totalWith1Component++;
        }
    }
    std::cout<<"RESULT"
            <<" method=cuckoo"
            <<" l="<<l
            <<" hfEvals="<<(double)hfEvals/(double)numIterations
            <<" tries="<<(double)totalTries / (double)numIterations
            <<" iterations="<<numIterations
            <<" spaceEstimate="<<log2((double)totalTries / (double)numIterations) / l + 1
            <<" pseudotrees="<< (double)totalPseudotrees/(double)numIterations
            <<" orientations="<< (double)totalOrientations/(double)numIterations
            <<" totalWith1Component="<< (double)totalWith1Component/(double)numIterations
            <<std::endl;
}

void testShockHashRot(size_t l) {
    assert(l < 64);
    bool actuallyRotate = true;
    size_t numIterations = l <= 30 ? 20000 : (l <= 43 ? 4000 : 1000);
    size_t totalTries = 0;
    size_t hfEvals = 0;
    for (size_t iteration = 0; iteration < numIterations; iteration++) {
        std::vector<uint64_t> keys;
        auto time = std::chrono::system_clock::now();
        long seed = std::chrono::duration_cast<std::chrono::milliseconds>(time.time_since_epoch()).count() * (iteration + 1);
        util::XorShift64 prng(seed);
        shockhash::TinyBinaryCuckooHashTable table(l);
        for (size_t i = 0; i < l; i++) {
            keys.push_back(prng());
            table.prepare(shockhash::HashedKey(keys.at(i)));
        }
        totalTries--;
        hfEvals -= 2 * l;
        uint64_t allSet = (1ul << l) - 1;
        bool found = false;
        mainLoop:while (!found) {
            totalTries++;
            hfEvals += 2 * l;

            shockhash::UnionFind unionFind(l);
            uint64_t a = 0;
            for (size_t i = 0; i < l/2; i++) {
                auto candidateCells = shockhash::TinyBinaryCuckooHashTable::getCandidateCells(shockhash::HashedKey(keys.at(i)), totalTries, l);
                a |= 1ull << candidateCells.cell1;
                a |= 1ull << candidateCells.cell2;
                if (!unionFind.unionIsStillPseudoforrest(candidateCells.cell1, candidateCells.cell2)) {
                    goto mainLoop; // First half does not work, so no reason to try second half
                }
            }
            std::vector<shockhash::TinyBinaryCuckooHashTable::CandidateCells> candidateCellsB(l);
            uint64_t b = 0;
            for (size_t i = l/2; i < l; i++) {
                auto candidateCells = shockhash::TinyBinaryCuckooHashTable::getCandidateCells(shockhash::HashedKey(keys.at(i)), totalTries, l);
                candidateCellsB.at(i) = candidateCells;
                b |= 1ull << candidateCells.cell1;
                b |= 1ull << candidateCells.cell2;
            }
            for (size_t r = 0; r < l; r++) {
                if ((a | rotate(l, b, r)) != allSet) {
                    continue;
                }
                if (!actuallyRotate && r != 0) {
                    continue;
                }
                auto unionFindCopy = unionFind;
                size_t i = l/2;
                for (; i < l; i++) {
                    auto candidateCells = candidateCellsB.at(i);
                    if (!unionFindCopy.unionIsStillPseudoforrest((candidateCells.cell1 + r) % l, (candidateCells.cell2 + r) % l)) {
                        break; // Try next rotation
                    }
                }
                if (i == l) {
                    // All were still pseudoforrests => Can be solved
                    found = true;
                    break;
                }
            }
        }
    }
    std::cout<<"RESULT"
             <<" method=cuckooRot"
             <<" l="<<l
             <<" hfEvals="<<(double)hfEvals/(double)numIterations
             <<" tries="<<(double)totalTries / (double)numIterations
             <<" iterations="<<numIterations
             <<" spaceEstimate="<<log2((double)totalTries * (double)(actuallyRotate ? l : 1) / (double)numIterations) / l + 1
             <<std::endl;
}

int main() {
    for (size_t l = 5; l <= 40; l++) {
        if (l <= 25) {
            testRotationFitting(l);
        }
        if (l <= 22) {
            testBruteForce(l);
        }
        testShockHash(l);
        testShockHashRot(l);
    }
}
