# ShockHash

Small, heavily overloaded cuckoo hash tables inside the [RecSplit](https://github.com/vigna/sux/blob/master/sux/function/RecSplit.hpp) framework.
Constructs very compact perfect hash functions significantly faster than previous approaches.

<img src="https://raw.githubusercontent.com/ByteHamster/ShockHash/main/plots.png" alt="Plots preview" />

### Library Usage

Clone (with submodules) this repo and add it to your `CMakeLists.txt`:

```
add_subdirectory(path/to/ShockHash)
target_link_libraries(YourTarget PRIVATE ShockHash)
```

### Licensing
ShockHash is licensed exactly like `libstdc++` (GPLv3 + GCC Runtime Library Exception), which essentially means you can use it everywhere, exactly like `libstdc++`.
You can find details in the [COPYING](/COPYING) and [COPYING.RUNTIME](/COPYING.RUNTIME) files.
