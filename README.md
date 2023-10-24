# ShockHash

Small, heavily overloaded cuckoo hash tables inside the [RecSplit](https://github.com/vigna/sux/blob/master/sux/function/RecSplit.hpp) framework.
Constructs very compact perfect hash functions significantly faster than previous approaches.

<img src="https://raw.githubusercontent.com/ByteHamster/ShockHash/main/plots.png" alt="Plots preview" />

### Library Usage

Clone (with submodules) this repo and add it to your `CMakeLists.txt`:

```cmake
add_subdirectory(path/to/ShockHash)
target_link_libraries(YourTarget PRIVATE ShockHash)
```

Then use one of the following classes:

- [ShockHash](https://github.com/ByteHamster/ShockHash/blob/main/include/ShockHash.h) is the original ShockHash algorithm integrated into the RecSplit framework.
- [SIMDShockHash](https://github.com/ByteHamster/ShockHash/blob/main/include/SIMDShockHash.hpp) is the SIMD-parallel version of the original ShockHash algorithm. Both ShockHash and the RecSplit framework are SIMD-parallelized. If this implementation is used on a machine without SIMD support, it is slower than the non-SIMD version because SIMD operations are emulated.
- [ShockHash2](https://github.com/ByteHamster/ShockHash/blob/main/include/ShockHash2.h) is the bipartite ShockHash algorithm. Only the inner ShockHash loop is SIMD-parallel, the RecSplit framework is not. If this implementation is used on a machine without SIMD support, the implementation uses sequential operations without explicitly emulating SIMD. To turn off SIMD, change to SIMD lanes of size 1 in [ShockHash2-internal.h](https://github.com/ByteHamster/ShockHash/blob/main/include/ShockHash2-internal.h).

Constructing a ShockHash perfect hash function is then straightforward:

```cpp
std::vector<std::string> keys = {"abc", "def", "123", "456"};
shockhash::ShockHash<30, false> shockHash(keys, 2000); // ShockHash base case size n=30, bucket size b=2000
std::cout << shockHash("abc") << " " << shockHash("def") << " "
          << shockHash("123") << " " << shockHash("456") << std::endl;
// Output: 1 3 2 0
```

We also give the base-case implementations without the RecSplit framework, which makes it easier to understand the main idea.

- Original [ShockHash](https://github.com/ByteHamster/ShockHash/blob/main/benchmark/bijections/ShockHash1.h).
- [Bipartite ShockHash](https://github.com/ByteHamster/ShockHash/blob/main/include/ShockHash2-internal.h). The outer loop that is also given in the pseudocode of the paper is given in `BijectionsShockHash2::findSeed`.

### Licensing
ShockHash is licensed exactly like `libstdc++` (GPLv3 + GCC Runtime Library Exception), which essentially means you can use it everywhere, exactly like `libstdc++`.
You can find details in the [COPYING](/COPYING) and [COPYING.RUNTIME](/COPYING.RUNTIME) files.

If you use ShockHash or bipartite ShockHash in an academic context or publication, please cite our papers:

```
@misc{lehmann2023shockhash,
      title={ShockHash: Towards Optimal-Space Minimal Perfect Hashing Beyond Brute-Force},
      author={Hans-Peter Lehmann and Peter Sanders and Stefan Walzer},
      year={2023},
      eprint={2308.09561},
      archivePrefix={arXiv},
      primaryClass={cs.DS}
}

@misc{lehmann2023bipartite,
      title={Bipartite ShockHash: Pruning ShockHash Search for Efficient Perfect Hashing},
      author={Hans-Peter Lehmann and Peter Sanders and Stefan Walzer},
      year={2023},
      eprint={2310.14959},
      archivePrefix={arXiv},
      primaryClass={cs.DS}
}
```
