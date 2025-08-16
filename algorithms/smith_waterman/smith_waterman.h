#ifndef SMITH_WATERMAN_H_
#define SMITH_WATERMAN_H_

#include <string>
#include <vector>

/// Holds the result of a Smith-Waterman alignment
struct SWResult {
  int score;                ///< Alignment score
  std::string aligned_seq1; ///< Aligned version of the first sequence
  std::string aligned_seq2; ///< Aligned version of the second sequence
};

/**
 * @brief Performs Smith-Waterman local sequence alignment.
 *
 * This function computes the optimal local alignment between two sequences
 * using the Smith-Waterman algorithm, returning the alignment score and
 * the aligned sequences.
 *
 * @param seq1 First input sequence
 * @param seq2 Second input sequence
 * @param match_score Score to assign for a match
 * @param mismatch_score Score to assign for a mismatch
 * @param gap_penalty Penalty for introducing a gap
 * @return SWResult containing the alignment score and the aligned sequences
 */
SWResult smith_waterman(const std::string &seq1, const std::string &seq2,
                        int match_score, int mismatch_score, int gap_penalty);

#endif // SMITH_WATERMAN_H_
