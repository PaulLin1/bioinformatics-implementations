#include "bioinformatics-pipeline/hamming_distance.h"
#include "bioinformatics-pipeline/cpu_hamming_distance.h"
// #include "bioinformatics-pipeline/cuda_hamming_distance.h"

#include <iostream>
#include <string>

int hamming_distance(const std::string& seq1, const std::string& seq2) {
    const int len1 = seq1.length();
    const int len2 = seq2.length();
    const int maxLen = std::max(len1, len2);

    // Pad sequences with '-'
    std::string padded1 = seq1 + std::string(maxLen - len1, '-');
    std::string padded2 = seq2 + std::string(maxLen - len2, '-');

    return cpu_hamming_distance(padded1, padded2);
}

// int main() {
//     std::string seq1 = "aTGACsd";
//     std::string seq2 = "ATGAC";

//     std::cout << "Hamming distance: " << hamming_distance(seq1, seq2) << std::endl;

//     return 0;
// }