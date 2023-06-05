cmake_minimum_required (VERSION 3.16)
cmake_policy(SET CMP0048 NEW)
project(ShockHash LANGUAGES CXX)

if(CMAKE_BUILD_TYPE STREQUAL "Release" AND CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_compile_options(-march=native)
endif()

# ---------------------------- Dependencies ----------------------------
if(NOT TARGET vectorclass)
    add_library(vectorclass INTERFACE)
    target_include_directories(vectorclass INTERFACE extlib/vectorclass)
endif()

if(NOT TARGET Sux)
    add_library(Sux INTERFACE)
    target_include_directories(Sux SYSTEM INTERFACE extlib/sux)
endif()

if(NOT TARGET tlx)
    set(TLX_INSTALL_INCLUDE_DIR tlx CACHE PATH "Workaround for TLX breaking the first cmake call")
    add_subdirectory(extlib/tlx)
endif()

if(NOT TARGET ips2ra)
    set(IPS2RA_DISABLE_PARALLEL ON CACHE PATH "ips2ra's FindTBB greps a file that does not exist in recent TBB versions")
    add_subdirectory(extlib/simpleRibbon)
    find_package(TBB)
    target_compile_options(ips2ra INTERFACE -D_REENTRANT)
    target_link_libraries(ips2ra INTERFACE pthread atomic TBB::tbb)
endif()

# ---------------------------- Library Setup ----------------------------

add_library(ShockHash INTERFACE)
target_include_directories(ShockHash INTERFACE include)
target_compile_features(ShockHash INTERFACE cxx_std_20)
target_link_libraries(ShockHash INTERFACE Sux SimpleRibbon)

include(${CMAKE_CURRENT_SOURCE_DIR}/extlib/cmake-findsse/FindSSE.cmake)
FindSSE()
if(SSE4_2_FOUND)
    add_library(ShockHashSIMD INTERFACE)
    target_link_libraries(ShockHashSIMD INTERFACE ShockHash vectorclass)
    target_compile_options(ShockHashSIMD INTERFACE -DSIMD -march=native)
endif()

# ---------------------------- Benchmarks ----------------------------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_library(BenchmarkUtils INTERFACE)
    target_include_directories(BenchmarkUtils INTERFACE benchmark)
    target_include_directories(BenchmarkUtils INTERFACE test)
    target_include_directories(BenchmarkUtils INTERFACE extlib/util/include)
    target_link_libraries(BenchmarkUtils INTERFACE tlx)

    # Warnings if this is the main project
    target_compile_options(ShockHash INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -frecord-gcc-switches>)

    add_executable(Benchmark benchmark/benchmark_construction.cpp)
    target_link_libraries(Benchmark PUBLIC BenchmarkUtils ShockHash)

    if(SSE4_2_FOUND)
        add_executable(BenchmarkSIMD benchmark/benchmark_construction.cpp)
        target_link_libraries(BenchmarkSIMD PUBLIC BenchmarkUtils ShockHashSIMD)
    endif()

    add_executable(NumHashEvals benchmark/numHashEvals.cpp)
    target_link_libraries(NumHashEvals PUBLIC BenchmarkUtils)
endif()
