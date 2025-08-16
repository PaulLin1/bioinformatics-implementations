#include "cpu_sw_create_scoring_matrix.h"

#include <algorithm>
#include <iostream>
#include <vector>

std::vector<int>
cpu_sw_create_scoring_matrix(const std::string &seq1, const std::string &seq2,
                             int match_score, int mismatch_score,
                             int gap_penalty, int rows, int cols) {
	std::vector<int> H(rows * cols, 0);

	for (int row = 1; row < rows; ++row) {
		for (int col = 1; col < cols; ++col) {
			int score =
			    (seq1[row - 1] == seq2[col - 1]) ? match_score : mismatch_score;

			int diag_score = H[(row - 1) * cols + (col - 1)] + score;
			int top_score = H[row * cols + (col - 1)] + gap_penalty;
			int left_score = H[(row - 1) * cols + col] + gap_penalty;

			H[row * cols + col] =
			    std::max({0, diag_score, top_score, left_score});
		}
	}

	return H;
}
