# Hamming Distance

---

## Overall Implementation Details

- If the 2 sequences are of different lengths, the extra letters are counted as mismatches
---

## CPU Implementation

Nothing to talk about. Really simple

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
