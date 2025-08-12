#include "smith_waterman.h"

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

std::tuple<std::vector<std::vector<int>>, std::pair<int, int>, int>
create_scoring_matrix(std::string seq1, std::string seq2, int match_score,
                      int mismatch_score, int gap_penalty) {
  // Create matrix of 0's that is length of seq1 + 1 by length of seq2 + 1
  int rows = seq1.length() + 1;
  int cols = seq2.length() + 1;
  std::vector<std::vector<int>> scoring_matrix(rows, std::vector<int>(cols, 0));

  std::pair<int, int> max_index = {0, 0};
  int max_value = 0;

  // Ignoring first col and row because those stay at 0
  for (int row = 1; row < scoring_matrix.size(); row++) {
    for (int col = 1; col < scoring_matrix[0].size(); col++) {
      int score =
          (seq1[row - 1] == seq2[col - 1]) ? match_score : mismatch_score;

      int diag_score = scoring_matrix[row - 1][col - 1] + score;
      int top_score = scoring_matrix[row][col - 1] + gap_penalty;
      int left_score = scoring_matrix[row - 1][col] + gap_penalty;

      scoring_matrix[row][col] =
          std::max({0, diag_score, top_score, left_score});

      // Track max value and its index
      if (scoring_matrix[row][col] > max_value) {
        max_value = scoring_matrix[row][col];
        max_index = {row, col};
      }
    }
  }

  return std::make_tuple(scoring_matrix, max_index, max_value);
}

std::pair<std::string, std::string>
traceback(const std::string &seq1, const std::string &seq2,
          const std::vector<std::vector<int>> &scoring_matrix,
          const std::pair<int, int> &best_index, int max_value, int match_score,
          int mismatch_score, int gap_penalty) {
  std::string seq1_sub;
  std::string seq2_sub;

  int current_score = max_value;
  std::pair<int, int> current_index = best_index;

  do {
    seq1_sub += seq1[current_index.first - 1];
    seq2_sub += seq2[current_index.second - 1];

    std::pair<int, int> diag_index = {current_index.first - 1,
                                      current_index.second - 1};
    std::pair<int, int> top_index = {current_index.first,
                                     current_index.second - 1};
    std::pair<int, int> left_index = {current_index.first - 1,
                                      current_index.second};

    int diag_score = scoring_matrix[diag_index.first][diag_index.second];
    int top_score = scoring_matrix[top_index.first][top_index.second];
    int left_score = scoring_matrix[left_index.first][left_index.second];

    if (current_score == diag_score + match_score ||
        current_score == diag_score + mismatch_score) {
      current_score = diag_score;
      current_index = diag_index;

    } else if (current_score == top_score + gap_penalty) {
      current_score = top_score;
      current_index = top_index;
    } else if (current_score == left_score + gap_penalty) {
      current_score = left_score;
      current_index = left_index;
    }
  } while (current_score != 0);

  // Reverse strings because be backtraced
  std::reverse(seq1_sub.begin(), seq1_sub.end());
  std::reverse(seq2_sub.begin(), seq2_sub.end());

  return std::make_pair(seq1_sub, seq2_sub);
}

SWResult smith_waterman(const std::string &seq1, const std::string &seq2,
                        int match_score, int mismatch_score, int gap_penalty) {
  auto [scoring_matrix, max_index, max_value] = create_scoring_matrix(
      seq1, seq2, match_score, mismatch_score, gap_penalty);
  auto [seq1_aligned, seq2_aligned] =
      traceback(seq1, seq2, scoring_matrix, max_index, max_value, match_score,
                mismatch_score, gap_penalty);

  return {max_value, seq1_aligned, seq2_aligned};
}
