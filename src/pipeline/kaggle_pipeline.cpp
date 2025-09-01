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
#include <cuda_runtime.h>

std::vector<std::string> chunk_sequence(const std::string& seq, int chunk_len, int overlap) {
    if (chunk_len <= 0 || overlap < 0 || overlap >= chunk_len)
        throw std::invalid_argument("Invalid chunk_len/overlap values");

    std::vector<std::string> chunks;
    int step = chunk_len - overlap;

    for (size_t i = 0; i < seq.size(); i += step) {
        int len = std::min(chunk_len, static_cast<int>(seq.size() - i));
        if (len > 0) chunks.push_back(seq.substr(i, len));
    }
    return chunks;
}

std::map<int, std::vector<std::string>> create_class_to_sequences(const std::string& file_name) {
    std::ifstream infile(file_name);
    if (!infile) throw std::runtime_error("Could not open database file: " + file_name);

    std::map<int, std::vector<std::string>> class_to_sequences;
    std::string line;
    std::getline(infile, line);  // skip header

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string seq, cls_str;
        if (std::getline(iss, seq, '\t') && std::getline(iss, cls_str)) {
            try {
                int cls = std::stoi(cls_str);
                class_to_sequences[cls].push_back(std::move(seq));
            } catch (...) {
                std::cerr << "Warning: skipping line with invalid class label: " << line << "\n";
            }
        }
    }
    return class_to_sequences;
}

std::map<int, float> detect_homology_gpu(
    const std::vector<std::string>& query_chunks,
    const std::vector<std::string>& db_chunks,
    const std::vector<int>& db_labels,
    const std::string& algorithm,
    int match_score, int mismatch_penalty, int gap_penalty) {

    std::map<int, float> class_to_score;
    std::map<int, int> class_counts;

    const int n_streams = 4;
    std::vector<cudaStream_t> streams(n_streams);
    for (int i = 0; i < n_streams; ++i) cudaStreamCreate(&streams[i]);

    for (size_t i = 0; i < db_chunks.size(); ++i) {
        int stream_id = i % n_streams;

        int result = 0;
        if (algorithm == "Smith-Waterman") {
            result = cuda_smith_waterman(query_chunks[0], db_chunks[i], match_score, mismatch_penalty, gap_penalty, streams[stream_id]).score;
        } else if (algorithm == "Needleman-Wunsch") {
            result = cuda_needleman_wunsch(query_chunks[0], db_chunks[i], match_score, mismatch_penalty, gap_penalty, streams[stream_id]).score;
        } else if (algorithm == "Hamming-Distance") {
            result = cuda_hamming_distance(query_chunks[0], db_chunks[i], streams[stream_id]);
        }

        class_to_score[db_labels[i]] += result;
        class_counts[db_labels[i]]++;
    }

    for (auto& s : streams) cudaStreamSynchronize(s);

    for (auto& [cls, score] : class_to_score) {
        score /= static_cast<float>(class_counts[cls]);
    }

    // Destroy streams
    for (auto& s : streams) cudaStreamDestroy(s);

    return class_to_score;
}

std::map<int, float> detect_homology(
    const std::string& query_sequence,
    const std::map<int, std::vector<std::string>>& class_to_sequences,
    const std::string& algorithm,
    int chunk_len,
    int overlap,
    bool use_cuda,
    int match_score = 2,
    int mismatch_penalty = -1,
    int gap_penalty = 2) {

    const auto query_chunks = chunk_sequence(query_sequence, chunk_len, overlap);

    if (!use_cuda) {
        // CPU path
        std::map<int, float> class_to_score;
        for (const auto& [label, sequences] : class_to_sequences) {
            std::vector<float> scores;
            for (const auto& seq : sequences) {
                auto seq_chunks = chunk_sequence(seq, chunk_len, overlap);
                for (const auto& q : query_chunks) {
                    for (const auto& s : seq_chunks) {
                        int result = 0;
                        if (algorithm == "Smith-Waterman") result = cpu_smith_waterman(q, s, match_score, mismatch_penalty, gap_penalty).score;
                        else if (algorithm == "Needleman-Wunsch") result = cpu_needleman_wunsch(q, s, match_score, mismatch_penalty, gap_penalty).score;
                        else if (algorithm == "Hamming-Distance") result = cpu_hamming_distance(q, s);
                        scores.push_back(result);
                    }
                }
            }
            class_to_score[label] = scores.empty() ? 0.0f : std::accumulate(scores.begin(), scores.end(), 0.0f) / scores.size();
        }
        return class_to_score;
    }

    std::vector<std::string> db_chunks_flat;
    std::vector<int> db_labels_flat;
    for (const auto& [label, sequences] : class_to_sequences) {
        for (const auto& seq : sequences) {
            auto chunks = chunk_sequence(seq, chunk_len, overlap);
            db_chunks_flat.insert(db_chunks_flat.end(), chunks.begin(), chunks.end());
            db_labels_flat.insert(db_labels_flat.end(), chunks.size(), label);
        }
    }

    return detect_homology_gpu(query_chunks, db_chunks_flat, db_labels_flat, algorithm, match_score, mismatch_penalty, gap_penalty);
}

struct Config {
    std::string db_file = "data/raw/kaggle_human.txt";
    std::string algorithm = "Hamming-Distance";
    int chunk_len = 64;
    int overlap = 16;
    bool use_cuda = true;
    std::string query_sequence = "ATGCGTACCTGAAGT";
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
        else if (key == "chunk_len") cfg.chunk_len = std::stoi(value);
        else if (key == "overlap") cfg.overlap = std::stoi(value);
        else if (key == "use_cuda") cfg.use_cuda = (value == "true" || value == "True" || value == "1");
        else if (key == "query_sequence") cfg.query_sequence = value;
        else throw std::runtime_error("Unknown argument: " + key);
    }

    return cfg;
}

int run_kaggle_pipeline(const Config &cfg) {
    auto class_to_sequences = create_class_to_sequences(cfg.db_file);

    auto scores = detect_homology(cfg.query_sequence,
                                  class_to_sequences,
                                  cfg.algorithm,
                                  cfg.chunk_len,
                                  cfg.overlap,
                                  cfg.use_cuda);

    auto best = std::max_element(
        scores.begin(), scores.end(),
        [](auto &a, auto &b) { return a.second < b.second; });

    std::cout << "Best match: class " << best->first
              << " with score " << best->second << "\n";
    return 0;
}

int main(int argc, char **argv) {
    try {
        Config cfg = parse_kwargs(argc, argv);
        return run_kaggle_pipeline(cfg);
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
