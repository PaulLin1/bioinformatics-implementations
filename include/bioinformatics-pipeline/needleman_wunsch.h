#ifndef NEEDLEMAN_WUNSCH_H_
#define NEEDLEMAN_WUNSCH_H_

#include <string>
#include <vector>

/// Holds the result of a Needleman-Wunsch alignment
struct NWResult {
	int score;                ///< Alignment score
	std::string aligned_seq1; ///< Aligned version of the first sequence
	std::string aligned_seq2; ///< Aligned version of the second sequence
};

/**
 * @brief Performs Needleman-Wunsch local sequence alignment.
 *
 * This function computes the optimal global alignment between two sequences
 * using the Needleman-Wunsch algorithm, returning the alignment score and
 * the aligned sequences.
 *
 * @param seq1 First input sequence
 * @param seq2 Second input sequence
 * @param match_score Score to assign for a match
 * @param mismatch_score Score to assign for a mismatch
 * @param gap_penalty Penalty for introducing a gap
 * @return NWResult containing the alignment score and the aligned sequences
 */
NWResult needleman_wunsch(const std::string &seq1, const std::string &seq2,
                        int match_score, int mismatch_score, int gap_penalty);

#endif // NEEDLEMAN_WUNSCH_H_
