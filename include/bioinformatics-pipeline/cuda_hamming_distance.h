#ifndef CUDA_HAMMING_DISTANCE_H_
#define CUDA_HAMMING_DISTANCE_H_

#include <string>
#include <cuda_runtime.h>

/**
 * @brief Finds Hamming Distance of 2 sequences.
 *
 * @param seq1 First input sequence
 * @param seq2 Second input sequence
 * @return Int containing the hamming distance
 */
int cuda_hamming_distance(const std::string& seq1, const std::string& seq2, cudaStream_t stream = 0);

#endif // CUDA_HAMMING_DISTANCE_H_
