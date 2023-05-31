#!/bin/bash
hostname
strings Benchmark | grep " -m"
strings BenchmarkSIMD | grep " -m"

# shellcheck disable=SC2086
for leafSize in $(seq 4 2 40); do
    for bucketSize in 5 50 500 2000; do
        params="--numObjects 500k --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"
        ./Benchmark $params
        ./BenchmarkSIMD $params
    done
done
