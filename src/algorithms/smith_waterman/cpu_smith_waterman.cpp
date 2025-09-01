#include "bioinformatics-pipeline/cpu_smith_waterman.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

std::vector<int>
cpu_sw_create_scoring_matrix(const std::string &seq1, const std::string &seq2,
                             int match_score, int mismatch_score,
                             int gap_penalty, int rows, int cols) {
    std::vector<int> H(rows * cols, 0);

    for (int row = 1; row < rows; ++row) {
        for (int col = 1; col < cols; ++col) {
            int score = (seq1[row - 1] == seq2[col - 1]) ? match_score : mismatch_score;

            int diag_score = H[(row - 1) * cols + (col - 1)] + score;
            int top_score  = H[row * cols + (col - 1)] + gap_penalty;
            int left_score = H[(row - 1) * cols + col] + gap_penalty;

            H[row * cols + col] = std::max({0, diag_score, top_score, left_score});
        }
    }

    return H;
}

static std::pair<std::string, std::string>
sw_traceback(const std::string &seq1, const std::string &seq2,
             const std::vector<int> &H, int rows, int cols,
             const std::pair<int, int> &best_index, int max_value,
             int match_score, int mismatch_score, int gap_penalty) {

    std::string seq1_sub, seq2_sub;
    int current_score = max_value;
    int i = best_index.first;
    int j = best_index.second;

    auto idx = [cols](int r, int c) { return r * cols + c; };

    while (current_score > 0) {
        if (i == 0 || j == 0) break;

        char a = seq1[i - 1];
        char b = seq2[j - 1];

        int diag_score = H[idx(i - 1, j - 1)];
        int top_score  = H[idx(i, j - 1)];
        int left_score = H[idx(i - 1, j)];

        int expected_diag_score = diag_score + ((a == b) ? match_score : mismatch_score);
        int expected_top_score  = top_score + gap_penalty;
        int expected_left_score = left_score + gap_penalty;

        if (current_score == expected_diag_score) {
            seq1_sub += a;
            seq2_sub += b;
            i -= 1; j -= 1;
            current_score = diag_score;
        }
        else if (current_score == expected_left_score) {
            seq1_sub += a;
            seq2_sub += '-';
            i -= 1;
            current_score = left_score;
        }
        else if (current_score == expected_top_score) {
            seq1_sub += '-';
            seq2_sub += b;
            j -= 1;
            current_score = top_score;
        }
        else {
            break; // safety
        }
    }

    std::reverse(seq1_sub.begin(), seq1_sub.end());
    std::reverse(seq2_sub.begin(), seq2_sub.end());

    return {seq1_sub, seq2_sub};
}

SWResult cpu_smith_waterman(const std::string &seq1, const std::string &seq2,
                        int match_score, int mismatch_score, int gap_penalty) {
    int rows = static_cast<int>(seq1.size()) + 1;
    int cols = static_cast<int>(seq2.size()) + 1;

    std::vector<int> H = cpu_sw_create_scoring_matrix(seq1, seq2,
                                                      match_score, mismatch_score,
                                                      gap_penalty, rows, cols);

    // Find max value and its position
    int max_val = 0, max_i = 0, max_j = 0;
    for (int i = 1; i < rows; ++i) {
        for (int j = 1; j < cols; ++j) {
            int v = H[i * cols + j];
            if (v > max_val) {
                max_val = v;
                max_i = i;
                max_j = j;
            }
        }
    }

    auto [aligned1, aligned2] = sw_traceback(seq1, seq2, H, rows, cols,
                                             {max_i, max_j}, max_val,
                                             match_score, mismatch_score, gap_penalty);
    return {max_val, aligned1, aligned2};

    // return {1, "a", "a"};
}

// int main() {
//     std::string seq1 = "ACACACTA";
//     std::string seq2 = "AGCACACA";

//     SWResult result = smith_waterman(seq1, seq2, 2, -1, -1);

//     std::cout << "Smith-Waterman result:\n";
//     std::cout << "Score: " << result.score << "\n";
//     std::cout << "Aligned seq1: " << result.aligned_seq1 << "\n";
//     std::cout << "Aligned seq2: " << result.aligned_seq2 << "\n";

//     return 0;
// }
