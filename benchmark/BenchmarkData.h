#pragma once
#include <vector>
#include <iostream>
#include <bytehamster/util/XorShift64.h>
#include <chrono>

std::vector<std::string> generateInputData(size_t N) {
    std::vector<std::string> inputData;
    inputData.reserve(N);
    auto time = std::chrono::system_clock::now();
    long seed = std::chrono::duration_cast<std::chrono::milliseconds>(time.time_since_epoch()).count();
    bytehamster::util::XorShift64 prng(seed);
    std::cout<<"Generating input data (Seed: "<<seed<<")"<<std::endl;
    char string[200];
    for (size_t i = 0; i < N; i++) {
        if ((i % (N/5)) == 0) {
            std::cout<<"\r\033[K"<<"Generating input: "<<100l*i/N<<"%"<<std::flush;
        }
        size_t length = 10 + prng((30 - 10) * 2);
        for (std::size_t k = 0; k < (length + sizeof(uint64_t))/sizeof(uint64_t); ++k) {
            ((uint64_t*) string)[k] = prng();
        }
        // Repair null bytes
        for (std::size_t k = 0; k < length; ++k) {
            if (string[k] == 0) {
                string[k] = 1 + prng(254);
            }
        }
        string[length] = 0;
        inputData.emplace_back(string, length);
    }
    std::cout<<"\r\033[K"<<"Input generation complete."<<std::endl;
    return inputData;
}
