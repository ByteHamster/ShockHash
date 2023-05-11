# HashFunc

Small, heavily overloaded cuckoo hash tables inside the [RecSplit](https://github.com/vigna/sux/blob/master/sux/function/RecSplit.hpp) framework.

### Library Usage

Clone (with submodules) this repo and add it to your `CMakeLists.txt`:

```
add_subdirectory(path/to/HashFunc)
target_link_libraries(YourTarget PRIVATE HashFunc)
```

### Licensing
HashFunc is licensed exactly like `libstdc++` (GPLv3 + GCC Runtime Library Exception), which essentially means you can use it everywhere, exactly like `libstdc++`.
You can find details in the [COPYING](/COPYING) and [COPYING.RUNTIME](/COPYING.RUNTIME) files.
