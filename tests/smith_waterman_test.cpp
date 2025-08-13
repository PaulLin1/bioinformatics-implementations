#include "smith_waterman.h"
#include "utils.h"
#include <gtest/gtest.h>

namespace {

constexpr int kMatchScore = 2;
constexpr int kMismatchScore = -1;
constexpr int kGapPenalty = -2;

TEST(SmithWatermanTest, CreateScoringMatrix) {
  const std::string seq1 = "GATTACA";
  const std::string seq2 = "GAATTACA";

  // Expected scoring matrix based on kMatchScore = 2, kMismatchScore = -1, kGapPenalty = -2
  // In case I switch params later
  std::vector<std::vector<int>> expected = {
      {  0,   0,   0,   0,   0,   0,   0,   0,   0 },
      {  0,   2,   0,   0,   0,   0,   0,   0,   0 },
      {  0,   0,   4,   2,   0,   0,   2,   0,   2 },
      {  0,   0,   2,   3,   4,   2,   0,   1,   0 },
      {  0,   0,   0,   1,   5,   6,   4,   2,   0 },
      {  0,   0,   2,   2,   3,   4,   8,   6,   4 },
      {  0,   0,   0,   1,   1,   2,   6,  10,   8 },
      {  0,   0,   2,   2,   0,   0,   4,   8,  12 }
  };

  auto [scoring_matrix, max_index, max_value] = create_scoring_matrix(seq1, seq2, kMatchScore, kMismatchScore, kGapPenalty);

  ASSERT_EQ(scoring_matrix.size(), expected.size()) << "Row count mismatch";

  for (size_t i = 0; i < expected.size(); ++i) {
    ASSERT_EQ(scoring_matrix[i].size(), expected[i].size()) << "Column count mismatch at row " << i;
    for (size_t j = 0; j < expected[i].size(); ++j) {
      EXPECT_EQ(scoring_matrix[i][j], expected[i][j])
          << "Mismatch at cell (" << i << "," << j << ")";
    }
  }

  EXPECT_EQ(max_index.first, 7);
  EXPECT_EQ(max_index.second, 8);
  EXPECT_EQ(max_value, 12);
}

// smith-waterman function tests
TEST(SmithWatermanTest, IdenticalSequences) {
  const std::string seq1 = "GATTACA";
  const std::string seq2 = "GATTACA";

  SWResult result =
      smith_waterman(seq1, seq2, kMatchScore, kMismatchScore, kGapPenalty);

  EXPECT_EQ(result.aligned_seq1, "GATTACA");
  EXPECT_EQ(result.aligned_seq2, "GATTACA");
  EXPECT_EQ(result.score, 7 * kMatchScore);
}

TEST(SmithWatermanTest, PartialMatch) {
  const std::string seq1 = "GATTACA";
  const std::string seq2 = "TTAC"; // Partial match in the middle

  SWResult result =
      smith_waterman(seq1, seq2, kMatchScore, kMismatchScore, kGapPenalty);

  EXPECT_EQ(result.aligned_seq1, "TTAC");
  EXPECT_EQ(result.aligned_seq2, "TTAC");
  EXPECT_EQ(result.score, 8);
}

TEST(SmithWatermanTest, SingleMismatch) {
  const std::string seq1 = "GATTACA";
  const std::string seq2 = "GACTACA";

  SWResult result =
      smith_waterman(seq1, seq2, kMatchScore, kMismatchScore, kGapPenalty);

  EXPECT_EQ(result.aligned_seq1, "GATTACA");
  EXPECT_EQ(result.aligned_seq2, "GACTACA");
  EXPECT_EQ(result.score, 6 * kMatchScore + kMismatchScore);
}

TEST(SmithWatermanTest, SingleGap) {
  const std::string seq1 = "GATTACA";
  const std::string seq2 = "GATACA";
  SWResult result = 
    smith_waterman(seq1, seq2, kMatchScore, kMismatchScore, kGapPenalty);
  EXPECT_EQ(result.aligned_seq1, "GATTACA");
  EXPECT_EQ(result.aligned_seq2, "GA-TACA");
  EXPECT_EQ(result.score, 6 * kMatchScore + kGapPenalty);
}

TEST(SmithWatermanTest, CompletelyDifferent) {
  const std::string seq1 = "AAAAAA";
  const std::string seq2 = "TTTTTT";

  SWResult result =
      smith_waterman(seq1, seq2, kMatchScore, kMismatchScore, kGapPenalty);

  EXPECT_EQ(result.aligned_seq1, "");
  EXPECT_EQ(result.aligned_seq2, "");
  EXPECT_EQ(result.score, 0);
}

TEST(SmithWatermanTest, OneEmptySequence) {
  const std::string seq1 = "";
  const std::string seq2 = "GATTACA";

  SWResult result =
      smith_waterman(seq1, seq2, kMatchScore, kMismatchScore, kGapPenalty);

  EXPECT_EQ(result.aligned_seq1, "");
  EXPECT_EQ(result.aligned_seq2, "");
  EXPECT_EQ(result.score, 0);
}

TEST(SmithWatermanTest, BothEmptySequences) {
  const std::string seq1 = "";
  const std::string seq2 = "";

  SWResult result =
      smith_waterman(seq1, seq2, kMatchScore, kMismatchScore, kGapPenalty);

  EXPECT_EQ(result.aligned_seq1, "");
  EXPECT_EQ(result.aligned_seq2, "");
  EXPECT_EQ(result.score, 0);
}

} // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
