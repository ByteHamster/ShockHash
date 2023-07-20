#!/bin/bash
hostname
strings recsplit_construction | grep " -m"
strings Benchmark | grep " -m"

# shellcheck disable=SC2086
for leafSize in $(seq 2 1 18); do
    for bucketSize in 100 500 1000 2000; do
        params="--numObjects 300k --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"
        ./recsplit_construction $params --leafMethod bruteforce
        ./recsplit_construction $params --leafMethod rotations
    done
done

# shellcheck disable=SC2086
for leafSize in $(seq 2 2 40); do
    for bucketSize in 100 500 1000 2000; do
        params="--numObjects 300k --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"
        ./Benchmark $params
        ./Benchmark $params --rotate
    done
done

# shellcheck disable=SC2086
for leafSize in $(seq 42 2 50); do
    for bucketSize in 100 500 1000 2000; do
        params="--numObjects 300k --numQueries 0 --leafSize $leafSize --bucketSize $bucketSize"
        ./Benchmark $params --rotate
    done
done
