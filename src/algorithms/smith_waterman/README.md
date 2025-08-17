# Smith Waterman

---

## Overall Implementation Details

- In my code, the scoring matrix is represented by the letter H, which is the standard variable used.
- My create_sw_scoring_matrix function specifies sw because needlman_wunsch uses it something similar.
- Everything matrix in row major.

---

## Naive CPU Implementation

Basic dynamic programming. Not much to talk about.

### Profiling Info with Perf on CPU

Regular Tests: 1.0s

1000 Genomes Project chromosome subsets (~2–4 GB FASTA files): too large to profile

NCBI RefSeq full organism collection: too large to profile

---

## Naive CUDA Implementation

The creation of the scoring matrix is parallelized.

### Profiling Info with Perf on H200

Regular Tests: 1.0s

1000 Genomes Project chromosome subsets (~2–4 GB FASTA files):

NCBI RefSeq full organism collection:
