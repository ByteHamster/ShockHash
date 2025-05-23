cmake_minimum_required (VERSION 3.25...4.0)
project(ShockHash
    DESCRIPTION "Various basic data structures"
    HOMEPAGE_URL "https://github.com/ByteHamster/ShockHash"
    VERSION 1.0
    LANGUAGES CXX)

if(TARGET ShockHash)
    return()
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if((CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo") AND PROJECT_IS_TOP_LEVEL)
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
    add_subdirectory(extlib/tlx SYSTEM EXCLUDE_FROM_ALL)
endif()

if(NOT TARGET simpleRibbon)
    set(IPS2RA_DISABLE_PARALLEL ON CACHE PATH "ips2ra's FindTBB greps a file that does not exist in recent TBB versions")
    add_subdirectory(extlib/simpleRibbon SYSTEM EXCLUDE_FROM_ALL)
    find_package(TBB)
    target_compile_options(ips2ra INTERFACE -D_REENTRANT)
    target_link_libraries(ips2ra INTERFACE pthread atomic TBB::tbb)
endif()

if(NOT TARGET Ips2raShockHashSorter)
    add_library(Ips2raShockHashSorter SHARED src/Sorter.cpp)
    target_compile_features(Ips2raShockHashSorter PRIVATE cxx_std_20)
    target_include_directories(Ips2raShockHashSorter PRIVATE include)
    target_link_libraries(Ips2raShockHashSorter PUBLIC ips2ra tlx Sux)
    target_compile_options(Ips2raShockHashSorter PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
endif()

if(NOT TARGET ByteHamsterUtil)
    add_subdirectory(extlib/util EXCLUDE_FROM_ALL)
endif()

# ---------------------------- Library Setup ----------------------------

add_library(ShockHash2Precompiled SHARED src/ShockHash2-precompiled.cpp)
target_compile_features(ShockHash2Precompiled PRIVATE cxx_std_20)
target_include_directories(ShockHash2Precompiled PRIVATE include)
target_link_libraries(ShockHash2Precompiled PUBLIC ips2ra tlx Sux ByteHamsterUtil)
target_compile_options(ShockHash2Precompiled PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)

add_library(ShockHash INTERFACE)
target_include_directories(ShockHash INTERFACE include)
target_compile_features(ShockHash INTERFACE cxx_std_20)
target_link_libraries(ShockHash INTERFACE Sux SimpleRibbon Ips2raShockHashSorter ByteHamsterUtil ShockHash2Precompiled)

include(${CMAKE_CURRENT_SOURCE_DIR}/extlib/cmake-findsse/FindSSE.cmake)
FindSSE()
if(SSE4_2_FOUND)
    add_library(ShockHashSIMD INTERFACE)
    target_link_libraries(ShockHashSIMD INTERFACE ShockHash vectorclass)
    target_compile_options(ShockHashSIMD INTERFACE -DSIMD -march=native)
endif()

# ---------------------------- Benchmarks ----------------------------
if(PROJECT_IS_TOP_LEVEL)
    add_library(BenchmarkUtils INTERFACE)
    target_include_directories(BenchmarkUtils INTERFACE benchmark)
    target_include_directories(BenchmarkUtils INTERFACE test)
    target_link_libraries(BenchmarkUtils INTERFACE tlx ByteHamster::Util)

    # Warnings if this is the main project
    target_compile_options(ShockHash INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Werror -Wno-error=stringop-overflow -frecord-gcc-switches>)

    add_executable(Benchmark benchmark/benchmark_construction.cpp)
    target_link_libraries(Benchmark PUBLIC BenchmarkUtils ShockHash)

    if(SSE4_2_FOUND)
        add_executable(BenchmarkSIMD benchmark/benchmark_construction.cpp)
        target_link_libraries(BenchmarkSIMD PUBLIC BenchmarkUtils ShockHashSIMD)
    endif()

    add_executable(NumHashEvals benchmark/numHashEvals.cpp)
    target_link_libraries(NumHashEvals PUBLIC BenchmarkUtils ShockHash)

    add_executable(Bijections benchmark/bijections.cpp)
    target_link_libraries(Bijections PUBLIC ShockHash BenchmarkUtils)

    add_executable(GolombMemoTuner benchmark/golombMemoTuner.cpp)
    target_link_libraries(GolombMemoTuner PUBLIC ShockHash BenchmarkUtils)
endif()
