/*
Pipeline for the kaggle dataset
Homology detection: Given a query sequence, align it against all sequences
in the database (class-labeled).
Data is .txt file with tab spacing
*/

#include "bioinformatics-pipeline/smith_waterman.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <map>
#include <algorithm>
#include <stdexcept>

// Split sequence into overlapping chunks (for GPU)
std::vector<std::string> chunk_sequence(const std::string& seq, int chunk_len, int overlap) {
    std::vector<std::string> chunks;

    if (chunk_len <= 0 || overlap < 0 || overlap >= chunk_len) {
        throw std::invalid_argument("Invalid chunk_len/overlap values");
    }

    int step = chunk_len - overlap;
    for (size_t i = 0; i < seq.size(); i += step) {
        int len = std::min(chunk_len, static_cast<int>(seq.size() - i));
        if (len > 0) {  // Avoid empty chunks
            chunks.push_back(seq.substr(i, len));
        }
    }

    return chunks;
}

// Load database: map from class to list of sequences
std::map<int, std::vector<std::string>> create_class_to_sequences(const std::string& file_name) {
    std::ifstream infile(file_name);
    if (!infile) {
        std::cerr << "Error: Could not open file " << file_name << "\n";
        return {};
    }
    
    std::map<int, std::vector<std::string>> class_to_sequences;
    std::string line;

    // Skip the header
    std::getline(infile, line);

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string sequence, cls_str;

        if (std::getline(iss, sequence, '\t') && std::getline(iss, cls_str)) {
            try {
                int cls = std::stoi(cls_str);
                class_to_sequences[cls].push_back(std::move(sequence));
            } catch (const std::invalid_argument&) {
                std::cerr << "Warning: Invalid class label, skipping line: " << line << "\n";
            }
        } else {
            std::cerr << "Warning: Skipping malformed line: " << line << "\n";
        }
    }

    return class_to_sequences;
}

// Compute similarity of query sequence against all database sequences
std::map<int, float> detect_homology(
    const std::string& query_sequence,
    const std::map<int, std::vector<std::string>>& class_to_sequences,
    const std::string& algorithm,
    int chunk_len,
    int overlap,
    bool CUDA,
    bool verbose = false,
    int match_score = 2,
    int mismatch_penalty = -1,
    int gap_penalty = 2) {
    
    // Chunk query sequence
    const auto query_chunks = chunk_sequence(query_sequence, chunk_len, overlap);

    if (query_chunks.empty()) {
        std::cerr << "Warning: Query sequence shorter than chunk length.\n";
    }

    std::map<int, float> class_to_score;

    for (const auto& [label, sequences] : class_to_sequences) {
        std::vector<float> all_scores;

        for (const auto& sequence : sequences) {
            const auto seq_chunks = chunk_sequence(sequence, chunk_len, overlap);

            if (seq_chunks.empty()) {
                if (verbose) std::cerr << "Warning: Sequence too short, skipping.\n";
                continue;
            }

            // Compare all query chunks to all sequence chunks
            for (const auto& q_chunk : query_chunks) {
                for (const auto& s_chunk : seq_chunks) {
                    SWResult result;
                    if (algorithm == "Smith-Waterman") {
                        result = smith_waterman(q_chunk, s_chunk,
                                                match_score, mismatch_penalty, gap_penalty);
                    } else {
                        throw std::invalid_argument("Invalid algorithm");
                    }
                    all_scores.push_back(static_cast<float>(result.score));
                }
            }
        }

        // Aggregate scores
        if (!all_scores.empty()) {
            float sum = 0.0f;
            for (float score : all_scores) sum += score;
            class_to_score[label] = sum / all_scores.size();
        } else {
            class_to_score[label] = 0.0f;
        }

        if (verbose) {
            std::cout << "Class " << label << " mean score: " << class_to_score[label] << "\n";
        }
    }

    return class_to_score;
}


int main() {
    const std::string db_file = "data/raw/kaggle_human.txt";
    const std::string algorithm = "Smith-Waterman";
    const int chunk_len = 64;
    const int overlap = 16;
    const bool CUDA = false;

    const std::string query_sequence = "ATGCGTACCTGAAGT";

    // Load database
    auto class_to_sequences = create_class_to_sequences(db_file);
    if (class_to_sequences.empty()) {
        std::cerr << "No database sequences loaded.\n";
        return 1;
    }

    // Detect homology
    auto class_scores = detect_homology(query_sequence, class_to_sequences,
                                        algorithm, chunk_len, overlap, CUDA);
    if (class_scores.empty()) {
        std::cerr << "Failed to compute homology scores\n";
        return 1;
    }

    // Find best-matching class
    auto best = std::max_element(class_scores.begin(), class_scores.end(),
                                 [](const auto& a, const auto& b) {
                                     return a.second < b.second;
                                 });

    std::cout << "\nBest match: class " << best->first
              << " with score " << best->second << "\n";

    return 0;
}
