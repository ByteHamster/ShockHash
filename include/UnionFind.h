#pragma once

namespace shockhash {
class UnionFind {
        std::vector<size_t> parents;
        std::vector<uint8_t> isTree_; // Avoid bit fiddling inside vector<bool>
        size_t n;
    public:
        explicit UnionFind(size_t n) : n(n) {
            parents.resize(n);
            isTree_.resize(n);
            clear();
        }

        void clear() {
            for (size_t i = 0;i < n; i++) {
                parents[i] = i;
                isTree_[i] = true;
            }
        }

        size_t findRepresentative(size_t x) {
            assert(x < n);
            while (parents[x] != x) {
                //size_t xOld = x;
                x = parents[x];
                //parents[xOld] = parents[x]; // Path compression is slower
            }
            return x;
        }

        bool unionIsStillPseudoforrest(size_t x, size_t y) {
            assert(x < n && y < n);
            size_t repr1 = findRepresentative(x);
            size_t repr2 = findRepresentative(y);
            parents[repr1] = repr2;
            //parents[x] = parents[y] = repr2; // Path compression is slower
            bool atLeastOneWasTree = isTree_[repr1] || isTree_[repr2];
            isTree_[repr2] = isTree_[repr1] && isTree_[repr2] && (repr1 != repr2);
            return atLeastOneWasTree;
        }
};
} // Namespace shockhash
