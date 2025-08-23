#include "bioinformatics-pipeline/cpu_hamming_distance.h"

#include <iostream>
#include <string>
#include <algorithm>

int cpu_hamming_distance(const std::string& seq1, const std::string& seq2) {
    int distance = 0;

    // Wrapper function makes them the same length
    for (size_t i = 0; i < seq1.length(); ++i) {
        char c1 = seq1[i];
        char c2 = seq2[i];
        if (c1 != c2) {
            distance++;
        }
    }

    return distance;
}


