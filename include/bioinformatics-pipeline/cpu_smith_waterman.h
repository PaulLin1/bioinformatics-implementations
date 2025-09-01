#ifndef CPU_SMITH_WATERMAN_H_
#define CPU_SMITH_WATERMAN_H_

#include <string>
#include <vector>

struct SWResult {
    int score;
    std::string aligned_seq1;
    std::string aligned_seq2;
};

/**
 * @brief Creates the Smith-Waterman scoring matrix on CPU.
 *
 * @param seq1 First input sequence
 * @param seq2 Second input sequence
 * @param match_score Score for a match
 * @param mismatch_score Score for a mismatch
 * @param gap_penalty Penalty for gaps
 * @param rows Number of rows for the scoring matrix (usually seq1.size() + 1)
 * @param cols Number of columns for the scoring matrix (usually seq2.size() +
 * 1)
 * @return 1D vector containing the scoring matrix in row-major order
 */

SWResult cpu_smith_waterman(const std::string &seq1, const std::string &seq2,
                            int match_score, int mismatch_score,
                            int gap_penalty);

#endif // CPU_SMITH_WATERMAN_H_
