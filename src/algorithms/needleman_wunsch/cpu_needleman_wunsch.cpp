#include "bioinformatics-pipeline/cpu_needleman_wunsch.h"

#include <vector>
#include <algorithm>
#include <string>
#include <iostream>

std::vector<int> cpu_nw_create_scoring_matrix(const std::string &seq1,
                                              const std::string &seq2,
                                              int match_score,
                                              int mismatch_score,
                                              int gap_penalty,
                                              int rows, int cols) {
    std::vector<int> H(rows * cols, 0);

    // Initialize gap penalties
    for (int i = 1; i < rows; ++i) H[i * cols + 0] = i * gap_penalty;
    for (int j = 1; j < cols; ++j) H[0 * cols + j] = j * gap_penalty;

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

static std::pair<std::string, std::string>
nw_traceback(const std::string &seq1, const std::string &seq2,
             const std::vector<int> &H, int rows, int cols,
             int match_score, int mismatch_score, int gap_penalty) {

    std::string seq1_sub, seq2_sub;
    int i = rows - 1;
    int j = cols - 1;

    auto idx = [cols](int r, int c) { return r * cols + c; };

    while (i > 0 || j > 0) {
        char a = (i > 0) ? seq1[i - 1] : '-';
        char b = (j > 0) ? seq2[j - 1] : '-';
        int score_current = H[idx(i, j)];

        if (i > 0 && j > 0) {
            int score_diag = H[idx(i - 1, j - 1)] + ((a == b) ? match_score : mismatch_score);
            if (score_current == score_diag) {
                seq1_sub += a;
                seq2_sub += b;
                --i; --j;
                continue;
            }
        }

        if (i > 0) {
            int score_up = H[idx(i - 1, j)] + gap_penalty;
            if (score_current == score_up) {
                seq1_sub += a;
                seq2_sub += '-';
                --i;
                continue;
            }
        }

        if (j > 0) {
            int score_left = H[idx(i, j - 1)] + gap_penalty;
            if (score_current == score_left) {
                seq1_sub += '-';
                seq2_sub += b;
                --j;
                continue;
            }
        }

        break; // safety fallback
    }

    std::reverse(seq1_sub.begin(), seq1_sub.end());
    std::reverse(seq2_sub.begin(), seq2_sub.end());
    return {seq1_sub, seq2_sub};
}

NWResult cpu_needleman_wunsch(const std::string &seq1, const std::string &seq2,
                          int match_score, int mismatch_score, int gap_penalty) {
    int rows = static_cast<int>(seq1.size()) + 1;
    int cols = static_cast<int>(seq2.size()) + 1;

    std::vector<int> H = cpu_nw_create_scoring_matrix(seq1, seq2,
                                                      match_score, mismatch_score,
                                                      gap_penalty, rows, cols);

    auto [aligned1, aligned2] = nw_traceback(seq1, seq2, H, rows, cols,
                                             match_score, mismatch_score, gap_penalty);

    int final_score = H[rows * cols - 1];
    return {final_score, aligned1, aligned2};
}

// // Example usage
// int main() {
//     std::string seq1 = "ACACACTA";
//     std::string seq2 = "AGCACACA";

//     NWResult result = needleman_wunsch(seq1, seq2, 2, -1, -1);

//     std::cout << "Score: " << result.score << "\n";
//     std::cout << "Aligned seq1: " << result.aligned_seq1 << "\n";
//     std::cout << "Aligned seq2: " << result.aligned_seq2 << "\n";

//     return 0;
// }
