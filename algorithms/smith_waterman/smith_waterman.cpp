
#include "smith_waterman.h"
#include "cpu_sw_create_scoring_matrix.h"
#include "cuda_sw_create_scoring_matrix.h"

#include <algorithm>
#include <iostream>
#include <utility>

static std::pair<std::string, std::string>
traceback(const std::string &seq1, const std::string &seq2,
          const std::vector<int> &H, int rows, int cols,
          const std::pair<int, int> &best_index, int max_value, int match_score,
          int mismatch_score, int gap_penalty) {

	std::string seq1_sub;
	std::string seq2_sub;

	int current_score = max_value;
	int i = best_index.first;
	int j = best_index.second;

	auto idx = [cols](int r, int c) { return r * cols + c; };

	while (current_score > 0) {
		// Indices for neighbors
		int diag_i = i - 1;
		int diag_j = j - 1;
		int top_i = i;
		int top_j = j - 1;
		int left_i = i - 1;
		int left_j = j;

		char a = seq1[i - 1];
		char b = seq2[j - 1];

		int diag_score =
		    (diag_i >= 0 && diag_j >= 0) ? H[idx(diag_i, diag_j)] : 0;
		int top_score = (top_j >= 0) ? H[idx(top_i, top_j)] : 0;
		int left_score = (left_i >= 0) ? H[idx(left_i, left_j)] : 0;

		int expected_diag_score =
		    diag_score + ((a == b) ? match_score : mismatch_score);
		int expected_top_score = top_score + gap_penalty;
		int expected_left_score = left_score + gap_penalty;

		if (current_score == expected_diag_score) {
			seq1_sub += a;
			seq2_sub += b;
			i = diag_i;
			j = diag_j;
			current_score = diag_score;

		} else if (current_score == expected_left_score) {
			seq1_sub += a;
			seq2_sub += '-';
			i = left_i;
			j = left_j;
			current_score = left_score;

		} else if (current_score == expected_top_score) {
			seq1_sub += '-';
			seq2_sub += b;
			i = top_i;
			j = top_j;
			current_score = top_score;

		} else {
			// No valid move found, break to avoid infinite loop
			break;
		}
	}

	// Reverse strings because the algorithm goes backwards
	std::reverse(seq1_sub.begin(), seq1_sub.end());
	std::reverse(seq2_sub.begin(), seq2_sub.end());

	return {seq1_sub, seq2_sub};
}

SWResult smith_waterman(const std::string &seq1, const std::string &seq2,
                        int match_score, int mismatch_score, int gap_penalty) {
	const int rows = static_cast<int>(seq1.size()) + 1;
	const int cols = static_cast<int>(seq2.size()) + 1;

	std::vector<int> H;

	#ifdef USE_CUDA
	H = cuda_sw_create_scoring_matrix(
		seq1, seq2, match_score, mismatch_score, gap_penalty, rows, cols);
	#else
	H = cpu_sw_create_scoring_matrix(
		seq1, seq2, match_score, mismatch_score, gap_penalty, rows, cols);
	#endif

	// Find the max value and matching index to perform traceback
	int max_val = 0;
	int max_i = 0, max_j = 0;
	for (int i = 1; i < rows; ++i) {
		for (int j = 1; j < cols; ++j) {
			// Row major format
			int v = H[i * cols + j];
			if (v > max_val) {
				max_val = v;
				max_i = i;
				max_j = j;
			}
		}
	}

	auto [aligned1, aligned2] =
	    traceback(seq1, seq2, H, rows, cols, {max_i, max_j}, max_val,
	              match_score, mismatch_score, gap_penalty);

	return {max_val, aligned1, aligned2};
}

// #include <iostream>

// int main() {
//     std::string seq1 = "ACACACTA";
//     std::string seq2 = "AGCACACA";

//     int match_score = 2;
//     int mismatch_score = -1;
//     int gap_penalty = -1;

//     SWResult result_cpu = smith_waterman(seq1, seq2, match_score,
//     mismatch_score, gap_penalty, false);

//     std::cout << "CPU Smith-Waterman result:\n";
//     std::cout << "Score: " << result_cpu.score << "\n";
//     std::cout << "Aligned seq1: " << result_cpu.aligned_seq1 << "\n";
//     std::cout << "Aligned seq2: " << result_cpu.aligned_seq2 << "\n\n";

//     SWResult result_cuda = smith_waterman(seq1, seq2, match_score,
//     mismatch_score, gap_penalty, true);

//     std::cout << "CUDA Smith-Waterman result:\n";
//     std::cout << "Score: " << result_cuda.score << "\n";
//     std::cout << "Aligned seq1: " << result_cuda.aligned_seq1 << "\n";
//     std::cout << "Aligned seq2: " << result_cuda.aligned_seq2 << "\n";

//     return 0;
// }