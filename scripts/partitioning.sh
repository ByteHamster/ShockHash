#!/bin/bash

hostname
for leafSize in $(seq 8 4 120); do
    params="--numObjects 100M --numQueries 100M --leafSize $leafSize"
    ./BenchmarkSIMD $params --shockhash2 --bucketSize 2000
    ./BenchmarkSIMD $params --shockhash2flat
done
