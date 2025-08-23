#include "bioinformatics-pipeline/needleman_wunsch.h"
#include "bioinformatics-pipeline/cpu_nw_create_scoring_matrix.h"
// #include "bioinformatics-pipeline/cuda_nw_create_scoring_matrix.h"

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>
#include <string>

static std::pair<std::string, std::string>
nw_traceback(const std::string &seq1, const std::string &seq2,
             const std::vector<int> &H, int rows, int cols,
             int match_score, int mismatch_score, int gap_penalty) {

    std::string seq1_sub;
    std::string seq2_sub;

    int i = rows - 1;
    int j = cols - 1;

    auto idx = [cols](int r, int c) { return r * cols + c; };

    while (i > 0 || j > 0) {
        char a = (i > 0) ? seq1[i - 1] : '-';
        char b = (j > 0) ? seq2[j - 1] : '-';

        int score_current = H[idx(i, j)];

        if (i > 0 && j > 0) {
            int score_diag = H[idx(i - 1, j - 1)] +
                             ((a == b) ? match_score : mismatch_score);
            if (score_current == score_diag) {
                seq1_sub += a;
                seq2_sub += b;
                --i;
                --j;
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

        break; // fallback safety
    }

    std::reverse(seq1_sub.begin(), seq1_sub.end());
    std::reverse(seq2_sub.begin(), seq2_sub.end());

    return {seq1_sub, seq2_sub};
}

NWResult needleman_wunsch(const std::string &seq1, const std::string &seq2,
                          int match_score, int mismatch_score, int gap_penalty) {
    const int rows = static_cast<int>(seq1.size()) + 1;
    const int cols = static_cast<int>(seq2.size()) + 1;

    std::vector<int> H = cpu_nw_create_scoring_matrix(
        seq1, seq2, match_score, mismatch_score, gap_penalty, rows, cols);

    auto [aligned1, aligned2] =
        nw_traceback(seq1, seq2, H, rows, cols, match_score, mismatch_score, gap_penalty);

    int final_score = H[rows * cols - 1]; // bottom-right cell

    return {final_score, aligned1, aligned2};
}

// // Example main()
// int main() {
//     std::string seq1 = "ACACACTA";
//     std::string seq2 = "AGCACACA";

//     int match_score = 2;
//     int mismatch_score = -1;
//     int gap_penalty = -1;

//     NWResult result = needleman_wunsch(seq1, seq2, match_score,
//                                        mismatch_score, gap_penalty);

//     std::cout << "Needleman-Wunsch result:\n";
//     std::cout << "Score: " << result.score << "\n";
//     std::cout << "Aligned seq1: " << result.aligned_seq1 << "\n";
//     std::cout << "Aligned seq2: " << result.aligned_seq2 << "\n";

//     return 0;
// }
