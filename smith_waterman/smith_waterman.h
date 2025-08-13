#ifndef SMITH_WATERMAN_H_
#define SMITH_WATERMAN_H_

#include <string>

// Holds the result of a Smith-Waterman alignment
struct SWResult {
  int score;
  std::string aligned_seq1;
  std::string aligned_seq2;
};

std::tuple<std::vector<std::vector<int>>, std::pair<int, int>, int>
create_scoring_matrix(std::string seq1, std::string seq2, int match_score,
                      int mismatch_score, int gap_penalty);


std::pair<std::string, std::string>
traceback(const std::string &seq1, const std::string &seq2,
          const std::vector<std::vector<int>> &scoring_matrix,
          const std::pair<int, int> &best_index, int max_value, int match_score,
          int mismatch_score, int gap_penalty);

SWResult smith_waterman(const std::string &seq1, const std::string &seq2,
                        int match_score, int mismatch_score, int gap_penalty);

#endif // SMITH_WATERMAN_H_