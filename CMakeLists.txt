cmake_minimum_required (VERSION 3.16)
cmake_policy(SET CMP0048 NEW)
project(ShockHash LANGUAGES CXX)

# ---------------------------- Dependencies ----------------------------
add_library(vectorclass INTERFACE)
target_include_directories(vectorclass INTERFACE extlib/vectorclass)

add_library(gcem INTERFACE)
target_include_directories(gcem INTERFACE extlib/gcem/include)

add_library(Sux INTERFACE)
target_include_directories(Sux SYSTEM INTERFACE extlib/sux)

set(IPS2RA_DISABLE_PARALLEL ON CACHE PATH "ips2ra's FindTBB greps a file that does not exist in recent TBB versions")
add_subdirectory(extlib/simpleRibbon)
find_package(TBB)
target_compile_options(ips2ra INTERFACE -D_REENTRANT)
target_link_libraries(ips2ra INTERFACE pthread atomic TBB::tbb)

# ---------------------------- Library Setup ----------------------------

add_library(ShockHash INTERFACE)
target_include_directories(ShockHash INTERFACE include)
target_compile_features(ShockHash INTERFACE cxx_std_20)
target_link_libraries(ShockHash INTERFACE Sux SimpleRibbon)

# ---------------------------- Benchmarks ----------------------------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	add_library(BenchmarkUtils INTERFACE)
	target_include_directories(BenchmarkUtils INTERFACE benchmark)
	target_include_directories(BenchmarkUtils INTERFACE test)
	target_include_directories(BenchmarkUtils INTERFACE extlib/util/include)

	set(TLX_INSTALL_INCLUDE_DIR tlx CACHE PATH "Workaround for TLX breaking the first cmake call")
	add_subdirectory(extlib/tlx EXCLUDE_FROM_ALL)
	target_link_libraries(BenchmarkUtils INTERFACE tlx)

	# Warnings if this is the main project
	target_compile_options(ShockHash INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -frecord-gcc-switches>)

	add_executable(Benchmark benchmark/benchmark_construction.cpp)
	target_link_libraries(Benchmark PUBLIC BenchmarkUtils ShockHash)

	add_executable(NumHashEvals benchmark/numHashEvals.cpp)
	target_link_libraries(NumHashEvals PUBLIC BenchmarkUtils)
endif()
