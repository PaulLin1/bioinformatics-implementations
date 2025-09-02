#include "bioinformatics-pipeline/cuda_smith_waterman.h"
#include <cuda_runtime.h>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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

void sw_init_buffers(SWBuffers& buf, int max_len1, int max_batch_size, int max_len2) {
    buf.max_len1 = max_len1;
    buf.max_batch_size = max_batch_size;
    buf.max_len2 = max_len2;

    CUDA_CHECK(cudaMalloc(&buf.d_seq1, max_len1 * sizeof(char)));
    CUDA_CHECK(cudaMalloc(&buf.d_seq2_padded, max_batch_size * max_len2 * sizeof(char)));
    CUDA_CHECK(cudaMalloc(&buf.d_seq2_lengths, max_batch_size * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&buf.d_H_batched, max_batch_size * (max_len1+1) * (max_len2+1) * sizeof(int)));

    buf.h_H_batched.resize(max_batch_size * (max_len1+1) * (max_len2+1));
}

void sw_free_buffers(SWBuffers& buf) {
    CUDA_CHECK(cudaFree(buf.d_seq1));
    CUDA_CHECK(cudaFree(buf.d_seq2_padded));
    CUDA_CHECK(cudaFree(buf.d_seq2_lengths));
    CUDA_CHECK(cudaFree(buf.d_H_batched));
}

static std::pair<std::string, std::string>
sw_traceback(const std::string &seq1, const std::string &seq2,
             const std::vector<int> &H, int rows, int cols,
             const std::pair<int, int> &best_index, int max_value,
             int match_score, int mismatch_score, int gap_penalty) {

    std::string seq1_sub, seq2_sub;
    int current_score = max_value;
    int i = best_index.first;
    int j = best_index.second;

    auto idx = [cols](int r, int c) { return r * cols + c; };

    while (current_score > 0) {
        if (i == 0 || j == 0) break;

        char a = seq1[i - 1];
        char b = seq2[j - 1];

        int diag_score = H[idx(i - 1, j - 1)];
        int top_score  = H[idx(i, j - 1)];
        int left_score = H[idx(i - 1, j)];

        int expected_diag_score = diag_score + ((a == b) ? match_score : mismatch_score);
        int expected_top_score  = top_score + gap_penalty;
        int expected_left_score = left_score + gap_penalty;

        if (current_score == expected_diag_score) {
            seq1_sub += a;
            seq2_sub += b;
            i -= 1; j -= 1;
            current_score = diag_score;
        }
        else if (current_score == expected_left_score) {
            seq1_sub += a;
            seq2_sub += '-';
            i -= 1;
            current_score = left_score;
        }
        else if (current_score == expected_top_score) {
            seq1_sub += '-';
            seq2_sub += b;
            j -= 1;
            current_score = top_score;
        }
        else {
            break; // safety
        }
    }

    std::reverse(seq1_sub.begin(), seq1_sub.end());
    std::reverse(seq2_sub.begin(), seq2_sub.end());

    return {seq1_sub, seq2_sub};
}

__global__ void sw_compute_diagonal_kernel_batched(
    const char* __restrict__ seq1,
    const char* __restrict__ seq2_padded,        
    const int*  __restrict__ seq2_lengths,  
    int* __restrict__ H_batched,            
    int len1, int max2, int rows, int cols,      
    int k, int i_start, int diag_len,
    int match_score, int mismatch_score, int gap_penalty)
{
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= diag_len) return;

    int b = blockIdx.y;              // which sequence in the batch
    const char* seq2 = seq2_padded + b * max2;
    int* H = H_batched + b * (rows * cols);
    int len2 = seq2_lengths[b];

    int i = i_start + t;
    int j = k - i;

    // guard (shouldn't happen if diag_len computed correctly)
    if (i <= 0 || j <= 0 || i >= rows || j >= cols) return;

    int a_idx = i - 1;          // index into seq1
    int b_idx = j - 1;          // index into seq2 (may be >= len2 → padding)

    int idx      = i * cols + j;
    int idx_diag = (i - 1) * cols + (j - 1);
    int idx_up   = (i - 1) * cols + j;
    int idx_left = i * cols + (j - 1);

    char a = seq1[a_idx];
    char bb = (b_idx < len2) ? seq2[b_idx] : 'N'; // padding as 'N'

    int score_sub = (a == bb) ? match_score : mismatch_score;

    int diag_score = H[idx_diag] + score_sub;
    int up_score   = H[idx_up]   + gap_penalty;
    int left_score = H[idx_left] + gap_penalty;

    int val = max(0, max(diag_score, max(up_score, left_score)));
    H[idx] = val;
}

void cuda_sw_create_scoring_matrix_batched(
    const std::string& seq1,
    const std::vector<std::string>& seq2_batch,
    int match_score, int mismatch_score, int gap_penalty,
    SWBuffers& buf)
{
    const int batch_size = (int)seq2_batch.size();
    const int len1 = (int)seq1.size();
    const int rows = len1 + 1;
    const int cols = buf.max_len2 + 1;

    if (len1 > buf.max_len1 || batch_size > buf.max_batch_size) {
        throw std::runtime_error("Batch size or seq length exceeds buffer limits");
    }

    std::vector<char> h_seq2_padded(batch_size * buf.max_len2, 'N');
    std::vector<int>  h_seq2_lengths(batch_size, 0);

    for (int b = 0; b < batch_size; ++b) {
        const auto& s = seq2_batch[b];
        int L = std::min((int)s.size(), buf.max_len2);
        std::copy_n(s.data(), L, &h_seq2_padded[b * buf.max_len2]);
        h_seq2_lengths[b] = L;
    }

    CUDA_CHECK(cudaMemcpy(buf.d_seq1, seq1.data(), len1 * sizeof(char), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(buf.d_seq2_padded, h_seq2_padded.data(),
                          h_seq2_padded.size() * sizeof(char), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(buf.d_seq2_lengths, h_seq2_lengths.data(),
                          batch_size * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(buf.d_H_batched, 0, batch_size * rows * cols * sizeof(int)));


    constexpr int TPB = 256; // threads per block

    for (int k = 2; k <= rows + cols - 2; ++k) {
        int i_start = max(1, k - (cols - 1));
        int i_end   = min(rows - 1, k - 1);
        int diag_len = i_end - i_start + 1;
        if (diag_len <= 0) continue;

        dim3 block(TPB);
        dim3 grid((diag_len + TPB - 1) / TPB, batch_size); // x = diag_len, y = batch_size

        sw_compute_diagonal_kernel_batched<<<grid, block>>>(
            buf.d_seq1,
            buf.d_seq2_padded,
            buf.d_seq2_lengths,
            buf.d_H_batched,
            len1, buf.max_len2, rows, cols,
            k, i_start, diag_len,
            match_score, mismatch_score, gap_penalty
        );

        CUDA_CHECK(cudaGetLastError());
    }

    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(buf.h_H_batched.data(), buf.d_H_batched,
                        batch_size * rows * cols * sizeof(int), cudaMemcpyDeviceToHost));
}

std::vector<SWResult> cuda_smith_waterman(
    const std::string& seq1,
    const std::vector<std::string>& seq2_batch,
    int match_score, int mismatch_score, int gap_penalty,
    SWBuffers& buf)
{
    const int MAX2 = buf.max_len2;
    const int batch_size = (int)seq2_batch.size();
    const int rows = (int)seq1.size() + 1;
    const int cols = MAX2 + 1;
    
    // Compute all scoring matrices on GPU using preallocated buffers
    cuda_sw_create_scoring_matrix_batched(
        seq1, seq2_batch, match_score, mismatch_score, gap_penalty, buf);


    const size_t elems = batch_size * rows * cols;
    buf.h_H_batched.resize(elems);
    cudaMemcpy(buf.h_H_batched.data(), buf.d_H_batched,
            elems * sizeof(int), cudaMemcpyDeviceToHost);

    // For each sequence, find max and traceback
    std::vector<SWResult> results;
    results.reserve(batch_size);

    for (int b = 0; b < batch_size; ++b) {
        const int base = b * rows * cols;
        const int* H = buf.h_H_batched.data() + base;

        int max_val = 0, max_i = 0, max_j = 0;
        for (int i = 1; i < rows; ++i) {
            const int row_off = i * cols;
            for (int j = 1; j < cols; ++j) {
                int v = H[row_off + j];
                if (v > max_val) {
                    max_val = v;
                    max_i = i;
                    max_j = j;
                }
            }
        }

        // Copy just this slice for traceback
        std::vector<int> H_slice(rows * cols);
        std::copy(buf.h_H_batched.begin() + base,
                buf.h_H_batched.begin() + base + rows * cols,
                H_slice.begin());

        auto [aligned1, aligned2] = sw_traceback(
            seq1, seq2_batch[b], H_slice, rows, cols,
            {max_i, max_j}, max_val,
            match_score, mismatch_score, gap_penalty);

        results.push_back({max_val, std::move(aligned1), std::move(aligned2)});
    }

    return results;
    // std::vector<SWResult> hi;
    // hi.push_back({1, "a", "a"});
    // return hi;
}


// int main() {
//     std::string seq1 = "ACACACTA";
//     std::string seq2 = "AGCACACA";
//     std::vector<std::string> hi;
//     hi.push_back(seq_2);

//     SWResult result = cuda_smith_waterman(seq1, hi, 2, -1, -1);

//     std::cout << "Smith-Waterman result:\n";
//     std::cout << "Score: " << result.score << "\n";
//     std::cout << "Aligned seq1: " << result.aligned_seq1 << "\n";
//     std::cout << "Aligned seq2: " << result.aligned_seq2 << "\n";

//     return 0;
// }