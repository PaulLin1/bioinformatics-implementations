/*
Tests work for both CPU and GPU
*/

#include "smith_waterman.h"

#include <gtest/gtest.h>

namespace {

constexpr int kMatchScore = 2;
constexpr int kMismatchScore = -1;
constexpr int kGapPenalty = -2;

// smith-waterman function tests
TEST(SmithWatermanTest, IdenticalSequences) {
	const std::string seq1 = "GATTACA";
	const std::string seq2 = "GATTACA";

	SWResult result = smith_waterman(seq1, seq2, kMatchScore, kMismatchScore,
	                                 kGapPenalty);

	EXPECT_EQ(result.aligned_seq1, "GATTACA");
	EXPECT_EQ(result.aligned_seq2, "GATTACA");
	EXPECT_EQ(result.score, 7 * kMatchScore);
}

TEST(SmithWatermanTest, PartialMatch) {
	const std::string seq1 = "GATTACA";
	const std::string seq2 = "TTAC"; // Partial match in the middle

	SWResult result = smith_waterman(seq1, seq2, kMatchScore, kMismatchScore,
	                                 kGapPenalty);

	EXPECT_EQ(result.aligned_seq1, "TTAC");
	EXPECT_EQ(result.aligned_seq2, "TTAC");
	EXPECT_EQ(result.score, 8);
}

TEST(SmithWatermanTest, SingleMismatch) {
	const std::string seq1 = "GATTACA";
	const std::string seq2 = "GACTACA";

	SWResult result = smith_waterman(seq1, seq2, kMatchScore, kMismatchScore,
	                                 kGapPenalty);

	EXPECT_EQ(result.aligned_seq1, "GATTACA");
	EXPECT_EQ(result.aligned_seq2, "GACTACA");
	EXPECT_EQ(result.score, 6 * kMatchScore + kMismatchScore);
}

TEST(SmithWatermanTest, SingleGap) {
	const std::string seq1 = "GATTACA";
	const std::string seq2 = "GATACA";
	SWResult result = smith_waterman(seq1, seq2, kMatchScore, kMismatchScore,
	                                 kGapPenalty);
	EXPECT_EQ(result.aligned_seq1, "GATTACA");
	EXPECT_EQ(result.aligned_seq2, "GA-TACA");
	EXPECT_EQ(result.score, 6 * kMatchScore + kGapPenalty);
}

TEST(SmithWatermanTest, CompletelyDifferent) {
	const std::string seq1 = "AAAAAA";
	const std::string seq2 = "TTTTTT";

	SWResult result = smith_waterman(seq1, seq2, kMatchScore, kMismatchScore,
	                                 kGapPenalty);

	EXPECT_EQ(result.aligned_seq1, "");
	EXPECT_EQ(result.aligned_seq2, "");
	EXPECT_EQ(result.score, 0);
}

TEST(SmithWatermanTest, OneEmptySequence) {
	const std::string seq1 = "";
	const std::string seq2 = "GATTACA";

	SWResult result = smith_waterman(seq1, seq2, kMatchScore, kMismatchScore,
	                                 kGapPenalty);

	EXPECT_EQ(result.aligned_seq1, "");
	EXPECT_EQ(result.aligned_seq2, "");
	EXPECT_EQ(result.score, 0);
}

TEST(SmithWatermanTest, BothEmptySequences) {
	const std::string seq1 = "";
	const std::string seq2 = "";

	SWResult result = smith_waterman(seq1, seq2, kMatchScore, kMismatchScore,
	                                 kGapPenalty);

	EXPECT_EQ(result.aligned_seq1, "");
	EXPECT_EQ(result.aligned_seq2, "");
	EXPECT_EQ(result.score, 0);
}

} // namespace

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
