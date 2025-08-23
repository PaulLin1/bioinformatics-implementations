#ifndef CPU_HAMMING_DISTANCE_H_
#define CPU_HAMMING_DISTANCE_H_

#include <string>

/**
 * @brief Finds Hamming Distance of 2 sequences.
 *
 * @param seq1 First input sequence
 * @param seq2 Second input sequence
 * @return Int containing the hamming distance
 */
int cpu_hamming_distance(const std::string& seq1, const std::string& seq2);

#endif // CPU_HAMMING_DISTANCE_H_
