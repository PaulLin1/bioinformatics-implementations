#include <iostream>
#include <string>
#include <algorithm>
#include <cuda_runtime.h>

/*
Macro for checking CUDA code
LLM recommended it. Might move to utils if it comes in handy.
*/
#ifndef CUDA_CHECK
#define CUDA_CHECK(expr)                                                       \
	do {                                                                       \
		cudaError_t _err = (expr);                                             \
		if (_err != cudaSuccess) {                                             \
			throw std::runtime_error(                                          \
			    std::string("CUDA error: ") + cudaGetErrorString(_err) +       \
			    " @ " + __FILE__ + ":" + std::to_string(__LINE__));            \
		}                                                                      \
	} while (0)
#endif

int cuda_hamming_distance(const std::string& seq1, const std::string& seq2) {
    const int len1 = seq1.length();
    const int len2 = seq2.length();
    const int maxLen = std::max(len1, len2);

    // Pad sequences with '-'
    std::string padded1 = seq1 + std::string(maxLen - len1, '-');
    std::string padded2 = seq2 + std::string(maxLen - len2, '-');

	// Pointers of sequences
	char *d_seq1 = nullptr, *d_seq2 = nullptr;
    int *d_result;
    int result = 0;

	// Allocate GPU mem
	CUDA_CHECK(cudaMalloc((void **)&d_seq1, max_len * sizeof(char)));
	CUDA_CHECK(cudaMalloc((void **)&d_seq2, max_len * sizeof(char)));
	CUDA_CHECK(cudaMalloc((void **)&d_res, sizeof(int)));

	// Copy seq1 and seq2 and init H to 0
	CUDA_CHECK(cudaMemcpy(d_seq1, padded1.data(), max_len * sizeof(char), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_seq2, padded2.data(), max_len * sizeof(char), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_seq2, &result, sizeof(int), cudaMemcpyHostToDevice));

    // Copy result back
    cudaMemcpy(&result, d_result, sizeof(int), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_seq1);
    cudaFree(d_seq2);
    cudaFree(d_result);

    return result;

}

// int main() {
//     std::string seq1 = "aTGACsd";
//     std::string seq2 = "ATGAC";

//     std::cout << "Hamming distance: " << cuda_hamming_distance(seq1, seq2) << std::endl;

//     return 0;
// }
