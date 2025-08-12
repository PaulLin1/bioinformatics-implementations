#include "smith_waterman.h"
#include <gtest/gtest.h>

namespace {

constexpr int kMatchScore = 2;
constexpr int kMismatchScore = -1;
constexpr int kGapPenalty = -2;

TEST(SmithWatermanTest, IdenticalSequences) {
  const std::string seq1 = "GATTACA";
  const std::string seq2 = "GATTACA";

  SWResult result =
      smith_waterman(seq1, seq2, kMatchScore, kMismatchScore, kGapPenalty);

  EXPECT_EQ(result.aligned_seq1, seq1);
  EXPECT_EQ(result.aligned_seq2, seq2);
  EXPECT_EQ(result.score,
            14);
}

} // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
