#include <vector>
#include <algorithm>

std::vector<int>
cpu_nw_create_scoring_matrix(const std::string &seq1, const std::string &seq2,
                             int match_score, int mismatch_score,
                             int gap_penalty, int rows, int cols) {
    std::vector<int> H(rows * cols, 0);

    // Initialization (gap penalties along top row and left col)
    for (int i = 1; i < rows; ++i) {
        H[i * cols + 0] = i * gap_penalty;
    }
    for (int j = 1; j < cols; ++j) {
        H[0 * cols + j] = j * gap_penalty;
    }

    // Fill scoring matrix
    for (int i = 1; i < rows; ++i) {
        for (int j = 1; j < cols; ++j) {
            int score = (seq1[i - 1] == seq2[j - 1]) ? match_score : mismatch_score;

            int diag_score = H[(i - 1) * cols + (j - 1)] + score;
            int top_score  = H[(i - 1) * cols + j] + gap_penalty;
            int left_score = H[i * cols + (j - 1)] + gap_penalty;

            H[i * cols + j] = std::max({diag_score, top_score, left_score});
        }
    }

    return H;
}