#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <utility>

// Smith-Waterman hyperparameters
constexpr int match =        2;
constexpr int mismatch =    -1;
constexpr int gap_penalty = -2;

void print_matrix(std::vector<std::vector<int> > matrix) {
    for (int row = 0; row < matrix.size(); row++) {
        for (int col = 0; col < matrix[0].size(); col++) {
            std::cout << matrix[row][col] << " ";
        }
        std::cout << std::endl;
    }
}

std::tuple<std::vector<std::vector<int>>, std::pair<int,int>, int> create_scoring_matrix(
    std::string sequence_1,
    std::string sequence_2) {
    // Create matrix of 0's that is length of sequence_1 + 1 by length of sequence_2 + 1
    int rows = sequence_1.length() + 1;
    int cols = sequence_2.length() + 1;
    std::vector<std::vector<int> > scoring_matrix(rows, std::vector<int>(cols, 0));

    std::pair<int, int> max_index = {0, 0};
    int max_value = 0;

    // Ignoring first col and row because those stay at 0
    for (int row = 1; row < scoring_matrix.size(); row++) {
        for (int col = 1; col < scoring_matrix[0].size(); col++) {
            int score = (sequence_1[row] == sequence_2[col]) ? match : mismatch;

            int top_left_score = scoring_matrix[row - 1][col - 1] + score;
            int top_score = scoring_matrix[row][col - 1] + gap_penalty;
            int left_score = scoring_matrix[row - 1][col] + gap_penalty;

            scoring_matrix[row][col] = std::max({0, top_left_score, top_score, left_score});

            // Track max value and its index
            if (scoring_matrix[row][col] > max_value) {
                max_value = scoring_matrix[row][col];
                max_index = {row, col};
            }
        }
    }

    return std::make_tuple(scoring_matrix, max_index, max_value);
}

std::pair<std::string, std::string> traceback(
    const std::string& sequence_1,
    const std::string& sequence_2,
    const std::vector<std::vector<int>>& scoring_matrix,
    const std::pair<int, int>& best_index,
    int max_value) {
    std::string sequence_1_subsqequence;
    std::string sequence_2_subsqequence;

    int current_score = max_value;
    std::pair<int, int> current_index = best_index;

    do {
        sequence_1_subsqequence += sequence_1[current_index.first - 1];
        sequence_2_subsqequence += sequence_1[current_index.second - 1];

        std::pair<int, int> top_left_index = {current_index.first - 1, current_index.second - 1};
        std::pair<int, int> top_index = {current_index.first, current_index.second - 1};
        std::pair<int, int> left_index = {current_index.first - 1, current_index.second};

        int top_left_score = scoring_matrix[top_left_index.first][top_left_index.second];
        int top_score = scoring_matrix[top_index.first][top_index.second];
        int left_score = scoring_matrix[left_index.first][left_index.second];

        if (current_score == top_left_score + match ||
            current_score == top_left_score + mismatch) {
            current_score = top_left_score;
            current_index = top_left_index;

        } else if (current_score == top_score + gap_penalty) {
            current_score = top_score;
            current_index = top_index;
        } else if (current_score == left_score + gap_penalty) {
            current_score = left_score;
            current_index = left_index;
            
        }
    } while (current_score != 0);

    // Reverse strings because be backtraced
    std::reverse(sequence_1_subsqequence.begin(), sequence_1_subsqequence.end());
    std::reverse(sequence_2_subsqequence.begin(), sequence_2_subsqequence.end());

    return std::make_pair(sequence_1_subsqequence, sequence_2_subsqequence);
}


int main() {
    std::string sequence_1 = "ACGTTGAC";
    std::string sequence_2 = "CGTTGA";

    auto [scoring_matrix, max_index, max_value] = create_scoring_matrix(sequence_1, sequence_2);
    auto [sequence_1_aligned, sequence_2_aligned] = traceback(sequence_1, sequence_2, scoring_matrix, max_index, max_value);

    std::cout << "Best local alignment score: " << max_value << std::endl;
    std::cout << "Sequence 1 Aligned: " << sequence_1_aligned << std::endl;
    std::cout << "Sequence 2 Aligned: " << sequence_2_aligned << std::endl;
}