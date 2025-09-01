#include "bioinformatics-pipeline/cuda_hamming_distance.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cuda_runtime.h>
#include <stdexcept>

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

__global__ void hamming_kernel(const char *seq1, const char *seq2, int *result, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        if (seq1[idx] != seq2[idx]) {
            atomicAdd(result, 1);
        }
    }
}

int cuda_hamming_distance(const std::string &seq1, const std::string &seq2, cudaStream_t stream) {
    const int len1 = seq1.length();
    const int len2 = seq2.length();
    const int maxLen = std::max(len1, len2);

    // Pad sequences with '-'
    std::string padded1 = seq1 + std::string(maxLen - len1, '-');
    std::string padded2 = seq2 + std::string(maxLen - len2, '-');

    // Device pointers
    char *d_seq1 = nullptr, *d_seq2 = nullptr;
    int *d_result = nullptr;
    int result = 0;

    // Allocate GPU memory asynchronously if stream != 0
    CUDA_CHECK(cudaMallocAsync((void **)&d_seq1, maxLen * sizeof(char), stream));
    CUDA_CHECK(cudaMallocAsync((void **)&d_seq2, maxLen * sizeof(char), stream));
    CUDA_CHECK(cudaMallocAsync((void **)&d_result, sizeof(int), stream));

    // Copy input sequences asynchronously
    CUDA_CHECK(cudaMemcpyAsync(d_seq1, padded1.data(), maxLen * sizeof(char), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_seq2, padded2.data(), maxLen * sizeof(char), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_result, &result, sizeof(int), cudaMemcpyHostToDevice, stream));

    // Launch kernel
    int threadsPerBlock = 256;
    int blocks = (maxLen + threadsPerBlock - 1) / threadsPerBlock;
    hamming_kernel<<<blocks, threadsPerBlock, 0, stream>>>(d_seq1, d_seq2, d_result, maxLen);
    CUDA_CHECK(cudaGetLastError());

    // Copy result back asynchronously
    CUDA_CHECK(cudaMemcpyAsync(&result, d_result, sizeof(int), cudaMemcpyDeviceToHost, stream));

    // Synchronize stream if non-zero
    if (stream != 0) CUDA_CHECK(cudaStreamSynchronize(stream));
    else CUDA_CHECK(cudaDeviceSynchronize());

    // Free device memory asynchronously if stream != 0
    CUDA_CHECK(cudaFreeAsync(d_seq1, stream));
    CUDA_CHECK(cudaFreeAsync(d_seq2, stream));
    CUDA_CHECK(cudaFreeAsync(d_result, stream));

    return result;
}
