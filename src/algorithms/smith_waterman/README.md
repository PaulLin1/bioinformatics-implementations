# Smith Waterman

---

## Overall Implementation Details

- In my code, the scoring matrix is represented by the letter H, which is the standard variable used.
- My create_sw_scoring_matrix function specifies sw because needlman_wunsch uses something similar. For clarity.
- Every matrix in row major.

---

## CPU Implementation

Basic dynamic programming. Not much to talk about.

### CPU Profiling Info with Perf

Regular Tests: 1.0s

Kaggle Dataset: 

NCBI RefSeq full organism collection: too large to profile

---

## CUDA Implementation

The creation of the scoring matrix is parallelized.

### H200 Profiling Info with Perf

Regular Tests: 1.0s

Kaggle Dataset: 

NCBI RefSeq full organism collection: too large to profile
