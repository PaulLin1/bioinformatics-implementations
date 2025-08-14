# Smith Waterman

More info about this algorithm can be found here: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

## Naive CPU Implementation
Basic dynamic programming. Not much to talk about.

### Profiling Info with Perf on CPU
Regular Tests: 1.0s
1000 Genomes Project chromosome subsets (~2–4 GB FASTA files): too large to profile
NCBI RefSeq full organism collection: too large to profile

## Naive CUDA Implementation
Basic dynamic programming. Not much to talk about. 

### Profiling Info with Perf on H200
Regular Tests: 1.0s
1000 Genomes Project chromosome subsets (~2–4 GB FASTA files): too large to profile
NCBI RefSeq full organism collection: too large to profile
