/*
Changed this a lot for the buffer stuff.
*/

#ifndef CUDA_SMITH_WATERMAN_H_
#define CUDA_SMITH_WATERMAN_H_

#include "bioinformatics-pipeline/cpu_smith_waterman.h"

#include <string>
#include <vector>
#include <cuda_runtime.h>

struct SWBuffers {
    char* d_seq1 = nullptr;
    char* d_seq2_padded = nullptr;
    int*  d_seq2_lengths = nullptr;
    int*  d_H_batched = nullptr;

    // Preallocated host-side buffer
    std::vector<int> h_H_batched;

    int max_len1 = 0;
    int max_batch_size = 0;
    int max_len2 = 0;
};


void sw_free_buffers(SWBuffers& buf);
void sw_init_buffers(SWBuffers& buf, int max_len1, int max_batch_size, int max_len2);

std::vector<SWResult> cuda_smith_waterman(
    const std::string& seq1,
    const std::vector<std::string>& seq2_batch,
    int match_score, int mismatch_score, int gap_penalty,
    SWBuffers& buf);

#endif // CUDA_SMITH_WATERMAN_H_
