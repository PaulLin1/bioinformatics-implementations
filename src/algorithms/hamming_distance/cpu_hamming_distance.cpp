#include "bioinformatics-pipeline/cpu_hamming_distance.h"

#include <string>
#include <algorithm>

int cpu_hamming_distance(const std::string& seq1, const std::string& seq2) {
    // Make sequences the same length by padding with '-'
    const size_t maxLen = std::max(seq1.length(), seq2.length());
    std::string padded1 = seq1 + std::string(maxLen - seq1.length(), '-');
    std::string padded2 = seq2 + std::string(maxLen - seq2.length(), '-');

    // Compute Hamming distance
    int distance = 0;
    for (size_t i = 0; i < maxLen; ++i) {
        if (padded1[i] != padded2[i]) {
            distance++;
        }
    }

    return distance;
}
