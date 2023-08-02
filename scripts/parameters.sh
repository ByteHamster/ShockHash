#!/bin/bash
hostname
strings Benchmark | grep " -m"
strings BenchmarkSIMD | grep " -m"

# shellcheck disable=SC2086
for leafSize in $(seq 4 2 40); do
    params="--numObjects 1M --numQueries 0 --leafSize $leafSize --bucketSize 2000"
    ./Benchmark $params
    ./BenchmarkSIMD $params
    ./BenchmarkSIMD $params --rotate
done
