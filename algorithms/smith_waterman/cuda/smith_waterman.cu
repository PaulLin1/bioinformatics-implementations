#include <cuda_runtime.h>
#include <iostream>

__global__ void diag_mul(const int* a, int* b, int n, int scalar) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_elements = n * n;
    if (idx >= total_elements) return;

    // Convert linear idx to anti-diagonal position
    int d = 0;
    int count = 0;
    while (true) {
        int elements_in_diag = (d < n) ? d + 1 : 2 * n - d - 1;
        if (idx < count + elements_in_diag) break;
        count += elements_in_diag;
        d++;
    }

    int offset_in_diag = idx - count;
    int i = (d < n) ? offset_in_diag : d - n + 1 + offset_in_diag;
    int j = d - i;

    b[i * n + j] = a[i * n + j] * scalar;
}

int main() {
    int n = 10, scalar = 5;
    int a[100], b[100];

    // Initialize 'a' matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i * n + j] = i; // row value
        }
    }

    int *d_a, *d_b;
    cudaMalloc(&d_a, n * n * sizeof(int));
    cudaMalloc(&d_b, n * n * sizeof(int));

    cudaMemcpy(d_a, a, n * n * sizeof(int), cudaMemcpyHostToDevice);

    int total_threads = n * n;
    int block_size = 128;
    int grid_size = (total_threads + block_size - 1) / block_size;

    diag_mul<<<grid_size, block_size>>>(d_a, d_b, n, scalar);
    cudaDeviceSynchronize();

    cudaMemcpy(b, d_b, n * n * sizeof(int), cudaMemcpyDeviceToHost);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << b[i * n + j] << " ";
        }
        std::cout << "\n";
    }

    cudaFree(d_a);
    cudaFree(d_b);
    return 0;
}
