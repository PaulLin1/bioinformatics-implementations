/*
Tests for Hamming Distance with padding
*/
#include "bioinformatics-pipeline/hamming_distance.h"

#include <gtest/gtest.h>
#include <string>

namespace {

// identical sequences
TEST(hamming_distanceTest, IdenticalSequences) {
    EXPECT_EQ(hamming_distance("GATTACA", "GATTACA"), 0);
}

// one mismatch in middle
TEST(hamming_distanceTest, SingleMismatch) {
    EXPECT_EQ(hamming_distance("GATTACA", "GACTACA"), 1);
}

// multiple mismatches, same length
TEST(hamming_distanceTest, MultipleMismatches) {
    EXPECT_EQ(hamming_distance("AAAAAA", "TTTTTT"), 6);
}

// different lengths, padding adds mismatches
TEST(hamming_distanceTest, DifferentLengthsExtraAtEnd) {
    EXPECT_EQ(hamming_distance("GATTACA", "GATTA"), 2);  
    // last two chars compared: 'C' vs '-', 'A' vs '-'
}

// different lengths, shorter first sequence
TEST(hamming_distanceTest, DifferentLengthsExtraInSecond) {
    EXPECT_EQ(hamming_distance("GAT", "GATTACA"), 4);  
    // "GAT---" vs "GATTACA" → mismatches at T vs T (match), then 3 vs 3 extra mismatches
}

// empty vs non-empty
TEST(hamming_distanceTest, OneEmptySequence) {
    EXPECT_EQ(hamming_distance("", "GATTACA"), 7);
    EXPECT_EQ(hamming_distance("GATTACA", ""), 7);
}

// both empty
TEST(hamming_distanceTest, BothEmptySequences) {
    EXPECT_EQ(hamming_distance("", ""), 0);
}

// completely different, same length
TEST(hamming_distanceTest, CompletelyDifferent) {
    EXPECT_EQ(hamming_distance("AAAA", "TTTT"), 4);
}

// mix of matches, mismatches, and padding
TEST(hamming_distanceTest, MixedCase) {
    EXPECT_EQ(hamming_distance("ACGT", "AGT"), 3);
    // Compare: A vs A (same), C vs G (diff), G vs T (diff), T vs - (diff)
    // total = 3
}

} // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
