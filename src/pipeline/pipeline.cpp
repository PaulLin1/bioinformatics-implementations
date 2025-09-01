#include "bioinformatics-pipeline/cpu_smith_waterman.h"
#include "bioinformatics-pipeline/cpu_hamming_distance.h"
#include "bioinformatics-pipeline/cpu_needleman_wunsch.h"

#include "bioinformatics-pipeline/cuda_smith_waterman.h"
#include "bioinformatics-pipeline/cuda_hamming_distance.h"
#include "bioinformatics-pipeline/cuda_needleman_wunsch.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <numeric>
#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <omp.h>
#include <iostream>
#include <chrono>

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

void run_pipeline_gpu(
    const std::string& query_sequence,
    const std::string& algorithm,
    const std::string& db_file,
    int match_score = 2,
    int mismatch_penalty = -1,
    int gap_penalty = 2) {
    std::ifstream infile(db_file);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << db_file << std::endl;
        return;
    }

    std::string line;
    std::string seq;
    seq.reserve(1 << 20);

    int sequences_processed = 0;

    std::vector<std::string> batch;
    const int BATCH_SIZE = 2048;

    auto start = std::chrono::high_resolution_clock::now();

    SWBuffers sw_buf;
    const int MAX_LEN1 = 500;      // maximum expected length of seq1
    const int MAX_BATCH = 2048;     // maximum batch size
    const int MAX_LEN2 = 500;      // maximum length of seq2 sequences

    sw_init_buffers(sw_buf, MAX_LEN1, MAX_BATCH, MAX_LEN2);


    while (getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!seq.empty()) {
                if (seq.size() <= 500) {
                    batch.push_back(seq);
                    sequences_processed++;
                }
                seq.clear();
            }
        } else {
            seq += line;
        }

        if (batch.size() == BATCH_SIZE) {
            // Process the whole batch at once
            std::vector<SWResult> hi = cuda_smith_waterman(query_sequence, batch,
                                    match_score, mismatch_penalty,
                                    gap_penalty, sw_buf);
            batch.clear();
        }

        if (sequences_processed >= 100000) break;
    }

    sw_free_buffers(sw_buf);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
}

void run_pipeline(
    const std::string& query_sequence,
    const std::string& algorithm,
    const std::string& db_file,
    int match_score = 2,
    int mismatch_penalty = -1,
    int gap_penalty = 2) {
    std::ifstream infile(db_file, std::ios::in);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << db_file << std::endl;
        return;
    }
    std::unique_ptr<char[]> bigbuf(new char[1 << 20]);
    infile.rdbuf()->pubsetbuf(bigbuf.get(), 1 << 20);

    std::string line;
    std::string seq;
    seq.reserve(1 << 20);

    int sequences_processed = 0;

    auto process_seq = [&](const std::string& current_seq) {
        if (current_seq.empty()) return;

        long long score = 0;

        if (algorithm == "Smith-Waterman") {
            int result = cpu_smith_waterman(
                query_sequence, current_seq,
                match_score, mismatch_penalty, gap_penalty
            ).score;
            score += result;
        } else if (algorithm == "Hamming-Distance") {
            int result = cpu_hamming_distance(
                query_sequence, current_seq
            );
            score += result;
        } else if (algorithm == "Needlman-Wunsch") {
            int result = cpu_needleman_wunsch(
                query_sequence, current_seq,
                match_score, mismatch_penalty, gap_penalty
            ).score;
            score += result;
        }

        std::cout << score << std::endl;
    };

    auto start = std::chrono::high_resolution_clock::now();

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (!line.empty() && line[0] == '>') {
            if (!seq.empty()) {
                process_seq(seq);
                seq.clear();
                sequences_processed++;
            }
        } else {
            seq.append(line);
        }

        if (sequences_processed > 100000) { break; }   
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
}

struct Config {
    std::string db_file = "/mnt/scratch/linpaul1/1uniprot_trembl.fasta";
    std::string algorithm = "Smith-Waterman";
    bool use_cuda = false;
    std::string query_sequence = "MKTFFVLLLFGVLTSASQAGDVEKNMKTFFVLLLFGMVLTSASQAGDVEKNLMAAHAGAVKAYTFFDLGQKGRTVQMGQMMKTFFVLLLMFGVLTSASQAGDVEKNLAAHAGAVKAYTFFDLGQKGRTVQGQLAAHAGAVKAYTFFDLGQKGRTVQGQMKTFFVLLLFGVLTSASQAGDVEKNLAAHAGAVKAYTFFDLGQMMKTFFVLLLFGVLTSASQAGDVEKNLAAHAGAVKAYTFFDLGQKGRTVQGQMKTFFVLLLFGVLTSASQAGDVEKNLAAHAGAVKAYTFFDLGQMKTFFVLLLFGVLTSASQAGDVEKNLAAHAGAVKAYTFFDLGQKGRTVQGQMKTFFVLLLFGVLTSASQAGDVEKNLAAHAGAVKAYTFFDLGQKGRTVQGQMKTFFVLLLFGVLTSASQAGDVEKNLAAHAGAVKAYTFFDLGQKGRTVQGQMKTFFVLLLFGVLTSASQAGDVEKNLAAHAGAVKAYTFFDLGQKGRTVQGQ";
};

Config parse_kwargs(int argc, char **argv) {
    Config cfg;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        auto eq_pos = arg.find('=');
        if (eq_pos == std::string::npos) {
            throw std::runtime_error("Invalid argument: " + arg);
        }
        std::string key = arg.substr(0, eq_pos);
        std::string value = arg.substr(eq_pos + 1);

        if (key == "db_file") cfg.db_file = value;
        else if (key == "algorithm") cfg.algorithm = value;
        else if (key == "use_cuda") cfg.use_cuda = (value == "true" || value == "True" || value == "1");
        else if (key == "query_sequence") cfg.query_sequence = value;
        else throw std::runtime_error("Unknown argument: " + key);
    }

    return cfg;
}

int main(int argc, char **argv) {
    try {
        Config cfg = parse_kwargs(argc, argv);
        if (!cfg.use_cuda) {
            run_pipeline(cfg.query_sequence,
                        cfg.algorithm,
                        cfg.db_file);
        } else {
            run_pipeline_gpu(cfg.query_sequence,
                             cfg.algorithm,
                             cfg.db_file);
        }

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
