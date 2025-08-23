/*
Tests for Needleman-Wunsch global alignment
*/
#include "bioinformatics-pipeline/needleman_wunsch.h"

#include <gtest/gtest.h>
#include <string>

namespace {

// identical sequences → no gaps, perfect score
TEST(NeedlemanWunschTest, IdenticalSequences) {
    NWResult result = needleman_wunsch("GATTACA", "GATTACA", 2, -1, -2);
    EXPECT_EQ(result.score, 7 * 2); // 7 matches * match_score(2)
    EXPECT_EQ(result.aligned_seq1, "GATTACA");
    EXPECT_EQ(result.aligned_seq2, "GATTACA");
}

// single mismatch
TEST(NeedlemanWunschTest, SingleMismatch) {
    NWResult result = needleman_wunsch("GATTACA", "GACTACA", 2, -1, -2);
    EXPECT_EQ(result.score, 6 * 2 + (-1)); // 6 matches, 1 mismatch
    EXPECT_EQ(result.aligned_seq1, "GATTACA");
    EXPECT_EQ(result.aligned_seq2, "GACTACA");
}

// completely different sequences
TEST(NeedlemanWunschTest, CompletelyDifferent) {
    NWResult result = needleman_wunsch("AAAA", "TTTT", 2, -1, -2);
    // best alignment = 4 mismatches, no benefit from gaps
    EXPECT_EQ(result.score, 4 * -1);
    EXPECT_EQ(result.aligned_seq1, "AAAA");
    EXPECT_EQ(result.aligned_seq2, "TTTT");
}

// one empty sequence
TEST(NeedlemanWunschTest, OneEmptySequence) {
    NWResult result = needleman_wunsch("GATTACA", "", 2, -1, -2);
    // score = 7 gaps * gap_penalty(-2) = -14
    EXPECT_EQ(result.score, -14);
    EXPECT_EQ(result.aligned_seq1, "GATTACA");
    EXPECT_EQ(result.aligned_seq2, "-------");

    result = needleman_wunsch("", "GATTACA", 2, -1, -2);
    EXPECT_EQ(result.score, -14);
    EXPECT_EQ(result.aligned_seq1, "-------");
    EXPECT_EQ(result.aligned_seq2, "GATTACA");
}

// both empty
TEST(NeedlemanWunschTest, BothEmpty) {
    NWResult result = needleman_wunsch("", "", 2, -1, -2);
    EXPECT_EQ(result.score, 0);
    EXPECT_EQ(result.aligned_seq1, "");
    EXPECT_EQ(result.aligned_seq2, "");
}

} // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
