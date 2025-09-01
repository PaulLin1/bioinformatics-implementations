#include "bioinformatics-pipeline/cuda_needleman_wunsch.h"

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cuda_runtime.h>

#ifndef CUDA_CHECK
#define CUDA_CHECK(expr)                                                       \
    do {                                                                       \
        cudaError_t _err = (expr);                                             \
        if (_err != cudaSuccess) {                                             \
            throw std::runtime_error(                                          \
                std::string("CUDA error: ") + cudaGetErrorString(_err) +       \
                " @ " + __FILE__ + ":" + std::to_string(__LINE__));            \
        }                                                                      \
    } while (0)
#endif

__global__ void nw_fill_diag(
    int* __restrict__ H,
    const char* __restrict__ s1,
    const char* __restrict__ s2,
    int rows, int cols,
    int match_score, int mismatch_score, int gap_penalty,
    int d
) {
    int i_start = max(1, d - (cols - 1));
    int i_end   = min(rows - 1, d - 1);
    int len = i_end - i_start + 1;
    if (len <= 0) return;

    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= len) return;

    int i = i_start + t;
    int j = d - i;

    int sub = (s1[i - 1] == s2[j - 1]) ? match_score : mismatch_score;

    int diag = H[(i - 1) * cols + (j - 1)] + sub;
    int up   = H[(i - 1) * cols + j]       + gap_penalty;
    int left = H[i * cols + (j - 1)]       + gap_penalty;

    H[i * cols + j] = max(diag, max(up, left));
}

std::vector<int>
cuda_nw_create_scoring_matrix(const std::string &seq1, const std::string &seq2,
                              int match_score, int mismatch_score,
                              int gap_penalty, int rows, int cols,
                              cudaStream_t stream) {
    std::vector<int> H_host(rows * cols, 0);
    for (int i = 1; i < rows; ++i) H_host[i * cols + 0] = i * gap_penalty;
    for (int j = 1; j < cols; ++j) H_host[0 * cols + j] = j * gap_penalty;

    int *d_H = nullptr;
    char *d_s1 = nullptr, *d_s2 = nullptr;

    CUDA_CHECK(cudaMallocAsync(&d_H,  rows * cols * sizeof(int), stream));
    CUDA_CHECK(cudaMallocAsync(&d_s1, (rows - 1) * sizeof(char), stream));
    CUDA_CHECK(cudaMallocAsync(&d_s2, (cols - 1) * sizeof(char), stream));

    CUDA_CHECK(cudaMemcpyAsync(d_H, H_host.data(),
                               rows * cols * sizeof(int),
                               cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_s1, seq1.data(),
                               (rows - 1) * sizeof(char),
                               cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_s2, seq2.data(),
                               (cols - 1) * sizeof(char),
                               cudaMemcpyHostToDevice, stream));

    const int TPB = 256;
    for (int d = 2; d <= rows + cols - 2; ++d) {
        int i_start = std::max(1, d - (cols - 1));
        int i_end   = std::min(rows - 1, d - 1);
        int len     = i_end - i_start + 1;
        if (len <= 0) continue;

        int blocks = (len + TPB - 1) / TPB;
        nw_fill_diag<<<blocks, TPB, 0, stream>>>(d_H, d_s1, d_s2,
                                                 rows, cols,
                                                 match_score, mismatch_score,
                                                 gap_penalty, d);
        CUDA_CHECK(cudaGetLastError());
    }

    // Copy back async
    CUDA_CHECK(cudaMemcpyAsync(H_host.data(), d_H,
                               rows * cols * sizeof(int),
                               cudaMemcpyDeviceToHost, stream));

    // Sync if non-default stream
    if (stream != 0) CUDA_CHECK(cudaStreamSynchronize(stream));
    else CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaFreeAsync(d_H, stream));
    CUDA_CHECK(cudaFreeAsync(d_s1, stream));
    CUDA_CHECK(cudaFreeAsync(d_s2, stream));

    return H_host;
}

static std::pair<std::string, std::string>
nw_traceback(const std::string &seq1, const std::string &seq2,
             const std::vector<int> &H, int rows, int cols,
             int match_score, int mismatch_score, int gap_penalty) {

    std::string seq1_sub, seq2_sub;
    int i = rows - 1;
    int j = cols - 1;
    auto idx = [cols](int r, int c) { return r * cols + c; };

    while (i > 0 || j > 0) {
        char a = (i > 0) ? seq1[i - 1] : '-';
        char b = (j > 0) ? seq2[j - 1] : '-';
        int score_current = H[idx(i, j)];

        if (i > 0 && j > 0) {
            int score_diag = H[idx(i - 1, j - 1)] + ((a == b) ? match_score : mismatch_score);
            if (score_current == score_diag) { seq1_sub += a; seq2_sub += b; --i; --j; continue; }
        }
        if (i > 0) { int score_up = H[idx(i - 1, j)] + gap_penalty;
                     if (score_current == score_up) { seq1_sub += a; seq2_sub += '-'; --i; continue; } }
        if (j > 0) { int score_left = H[idx(i, j - 1)] + gap_penalty;
                     if (score_current == score_left) { seq1_sub += '-'; seq2_sub += b; --j; continue; } }

        break;
    }

    std::reverse(seq1_sub.begin(), seq1_sub.end());
    std::reverse(seq2_sub.begin(), seq2_sub.end());
    return {seq1_sub, seq2_sub};
}

NWResult cuda_needleman_wunsch(const std::string &seq1, const std::string &seq2,
                               int match_score, int mismatch_score,
                               int gap_penalty, cudaStream_t stream) {
    int rows = static_cast<int>(seq1.size()) + 1;
    int cols = static_cast<int>(seq2.size()) + 1;

    std::vector<int> H = cuda_nw_create_scoring_matrix(seq1, seq2,
                                                       match_score, mismatch_score,
                                                       gap_penalty, rows, cols,
                                                       stream);

    auto [aligned1, aligned2] = nw_traceback(seq1, seq2, H, rows, cols,
                                             match_score, mismatch_score, gap_penalty);
    int final_score = H[rows * cols - 1];
    return {final_score, aligned1, aligned2};
}
