#ifndef SMITH_WATERMAN_H_
#define SMITH_WATERMAN_H_

#include <string>

// Holds the result of a Smith-Waterman alignment
struct SWResult {
  int score;
  std::string aligned_seq1;
  std::string aligned_seq2;
};

SWResult smith_waterman(const std::string &seq1, const std::string &seq2,
                        int match_score, int mismatch_score, int gap_penalty);

#endif // SMITH_WATERMAN_H_