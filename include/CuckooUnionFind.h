#pragma once
#include <cstdint>
#include <vector>

class CuckooUnionFind {
    private:
        static constexpr uint32_t EMPTY = ~uint32_t(0);
        std::vector<uint32_t> cells;
    public:
        explicit CuckooUnionFind(size_t numEntries) {
            cells.resize(numEntries);
        }

        void clear() {
            std::fill(cells.begin(), cells.end(), EMPTY);
        }

        bool unionIsStillPseudoforest(uint32_t cell1, uint32_t cell2) {
            uint32_t entry = cell1 ^ cell2;
            uint32_t origEntry = entry;
            if (cells[cell1] == EMPTY) {
                cells[cell1] = entry;
                return true;
            }
            if (cells[cell2] == EMPTY) {
                cells[cell2] = entry;
                return true;
            }
            uint32_t currentCell = cell1;
            uint32_t origCell = currentCell;

            do {
                uint32_t alternativeCell = entry ^ currentCell;
                std::swap(entry, cells[alternativeCell]);
                if (entry == EMPTY) {
                    return true;
                }
                currentCell = alternativeCell;
            } while (currentCell != origCell || entry != origEntry);
            return false;
        }
};
