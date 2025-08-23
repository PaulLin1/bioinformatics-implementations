#include "bioinformatics-pipeline/cuda_sw_create_scoring_matrix.h"

#include <algorithm>
#include <cuda_runtime.h>
#include <iostream>
#include <stdexcept>

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

// Cuda code
__global__ void sw_compute_diagonal_kernel(
    const char *__restrict__ seq1, const char *__restrict__ seq2,
    int *__restrict__ H, int rows, int cols, int k, int i_start, int diag_len,
    int match_score, int mismatch_score, int gap_penalty) {
	int t = blockIdx.x * blockDim.x + threadIdx.x;
	if (t >= diag_len)
		return;

	int i = i_start + t;
	int j = k - i;

	int a_idx = i - 1;
	int b_idx = j - 1;

	int idx = i * cols + j;
	int idx_diag = (i - 1) * cols + (j - 1);
	int idx_up = (i - 1) * cols + j;
	int idx_left = i * cols + (j - 1);

	int score_sub = (seq1[a_idx] == seq2[b_idx]) ? match_score : mismatch_score;

	int diag_score = H[idx_diag] + score_sub;
	int up_score = H[idx_up] + gap_penalty;
	int left_score = H[idx_left] + gap_penalty;

	int val = diag_score;
	if (up_score > val)
		val = up_score;
	if (left_score > val)
		val = left_score;
	if (val < 0)
		val = 0;

	H[idx] = val;
}

// Creates scoring matrix 1D
std::vector<int>
cuda_sw_create_scoring_matrix(const std::string &seq1, const std::string &seq2,
                              int match_score, int mismatch_score,
                              int gap_penalty, int rows, int cols) {
	/*
	Idk what to do rn for len and row/col.
	Both are needed.
	Should i pass in len or rows?
	Rows is used in the main sw func so I pass
	that in and recalc here for len1 and len2.
	Might switch
	*/
	const int len1 = rows - 1;
	const int len2 = rows - 1;

	// Pointers of sequences
	char *d_seq1 = nullptr, *d_seq2 = nullptr;
	// Pointer for the scoring matrix
	int *d_H = nullptr;

	// Allocate GPU mem
	CUDA_CHECK(cudaMalloc((void **)&d_seq1, len1 * sizeof(char)));
	CUDA_CHECK(cudaMalloc((void **)&d_seq2, len2 * sizeof(char)));
	CUDA_CHECK(cudaMalloc((void **)&d_H, rows * cols * sizeof(int)));

	// Copy seq1 and seq2 and init H to 0
	CUDA_CHECK(cudaMemcpy(d_seq1, seq1.data(), len1 * sizeof(char), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_seq2, seq2.data(), len2 * sizeof(char), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemset(d_H, 0, rows * cols * sizeof(int)));

	constexpr int TPB = 256; // Not optimized. Just random for now
	// Start looping from second diagonal
	// k is the diag number
	for (int k = 2; k <= rows + cols - 2; ++k) {
		int i_start = max(1, k - (cols - 1));
		int i_end = min(rows - 1, k - 1);
		int diag_len = i_end - i_start + 1;
		if (diag_len <= 0)
			continue;

		int blocks = (diag_len + TPB - 1) / TPB;
		// Launch the actual kernel
		sw_compute_diagonal_kernel<<<blocks, TPB>>>(
		    d_seq1, d_seq2, d_H, rows, cols, k, i_start, diag_len, match_score,
		    mismatch_score, gap_penalty);
		CUDA_CHECK(cudaGetLastError());
		CUDA_CHECK(cudaDeviceSynchronize());
	}

	// Copy everything back into a 1D matrix
	std::vector<int> H(rows * cols);
	CUDA_CHECK(cudaMemcpy(H.data(), d_H, rows * cols * sizeof(int),
	                      cudaMemcpyDeviceToHost));

	cudaFree(d_seq1);
	cudaFree(d_seq2);
	cudaFree(d_H);

	return H;
}
