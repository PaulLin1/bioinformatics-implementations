#ifndef CUDA_NEEDLEMAN_WUNSCH_H_
#define CUDA_NEEDLEMAN_WUNSCH_H_

#include "bioinformatics-pipeline/cpu_needleman_wunsch.h"

#include <string>
#include <vector>
#include <cuda_runtime.h>

/**
 * @brief Creates the Needleman-Wunsch scoring matrix on GPU.
 *
 * @param seq1 First input sequence
 * @param seq2 Second input sequence
 * @param match_score Score for a match
 * @param mismatch_score Score for a mismatch
 * @param gap_penalty Penalty for gaps
 * @param rows Number of rows for the scoring matrix (usually seq1.size() + 1)
 * @param cols Number of columns for the scoring matrix (usually seq2.size() +
 * 1)
 * @return 1D vector containing the scoring matrix in row-major order
 */
NWResult cuda_needleman_wunsch(const std::string &seq1, const std::string &seq2,
                             int match_score, int mismatch_score,
                             int gap_penalty, cudaStream_t stream = 0);

                            
#endif // CUDA_NEEDLEMAN_WUNSCH_H_