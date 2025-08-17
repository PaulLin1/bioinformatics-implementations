/*
Pipeline for the kaggle dataset
Class prediction using the algorithm scores
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

std::vector<std::string> chunk_sequence(const std::string& seq, int chunk_len, int overlap) {
    std::vector<std::string> chunks;

    if (chunk_len <= 0 || overlap < 0 || overlap >= chunk_len) {
        throw std::invalid_argument("Invalid chunk_len/overlap values");
    }

    for (size_t i = 0; i < seq.size(); i += (chunk_len - overlap)) {
        int len = std::min(chunk_len, static_cast<int>(seq.size() - i));
        chunks.push_back(seq.substr(i, len));
    }

    return chunks;
}

std::pair<int, std::string> pick_target_sequence(std::map<int, std::vector<std::string>>& data_by_class) {
    if (data_by_class.empty()) {
        throw std::runtime_error("No data to pick from!");
    }

    // Pick random key
    int key_index = std::rand() % data_by_class.size();
    auto it = data_by_class.begin();
    std::advance(it, key_index);

    auto& vec = it->second;
    if (vec.empty()) {
        throw std::runtime_error("Selected class has no sequences!");
    }

    // Pick random element from the vector
    int val_index = std::rand() % vec.size();
    std::string sequence = vec[val_index];

    // Remove the element from the vector
    vec.erase(vec.begin() + val_index);

    return {it->first, sequence};
}

std::map<int, std::vector<std::string>> create_class_to_sequences(std::string file_name) {
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
                class_to_sequences[cls].push_back(sequence);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Warning: Invalid class label, skipping line: " << line << "\n";
            }
        } else {
            std::cerr << "Warning: Skipping malformed line: " << line << "\n";
        }
    }

    return class_to_sequences;
}

void kaggle_pipeline(
    const std::string& file_name,
    const std::string& algorithm,
    int chunk_len,
    int overlap,
    bool CUDA,
    int match_score = 2,
    int mismatch_penalty =-1,
    int gap_penalty = 2) {
    /*
    Read in data and put them into a hashmap with the
    key being the class and the value being a list of the sequences
    in the class
    */
    std::map<int, std::vector<std::string>> class_to_sequences = create_class_to_sequences(file_name);
    if (class_to_sequences.empty()) {
        std::cerr << "No data loaded.\n";
    }

    // Get random value to compare against
    std::srand(0);
    // std::srand(std::time(nullptr)); // seed
    auto [target_sequence_class, target_sequence] = pick_target_sequence(class_to_sequences);
    std::cout << "Random value: " << target_sequence << " from class " << target_sequence_class << "\n";

    // Chunk target_sequence
    std::vector<std::string> target_seq_chunks = chunk_sequence(target_sequence, chunk_len, overlap);

    // Assigns a score to each class
    std::map<int, float> class_to_score;

    for (const auto& [label, sequences] : class_to_sequences) {
        int score_sum = 0;
        int class_length = sequences.size();

        // Preallocate chunk_scores vector once (max possible size)
        std::vector<int> chunk_scores;
        chunk_scores.reserve(5);

        // Temporary vector to hold chunks (reuse for each sequence)
        std::vector<std::string> current_seq_chunks;

        for (const auto& sequence : sequences) {
            // Reuse vector instead of creating new one
            current_seq_chunks.clear();
            current_seq_chunks = chunk_sequence(sequence, chunk_len, overlap);

            int num_of_iterations = std::min(current_seq_chunks.size(), target_seq_chunks.size());

            // Reuse chunk_scores vector
            chunk_scores.clear();

            if (algorithm == "Smith-Waterman") {
                for (int i = 0; i < num_of_iterations; i++) {
                    SWResult result_cpu = smith_waterman(target_seq_chunks[i], current_seq_chunks[i],
                                                         match_score, mismatch_penalty, gap_penalty);
                    chunk_scores.push_back(result_cpu.score);
                }
            } else {
                throw std::invalid_argument("Invalid algorithm");
            }

            // Add the best chunk score for this sequence
            score_sum += *std::max_element(chunk_scores.begin(), chunk_scores.end());
        }

        float average_score = static_cast<float>(score_sum) / class_length;
        class_to_score[label] = average_score;
        std::cout << "Class " << label << ": " << average_score << std::endl;
    }

    // Returns the class with the highest avg score
    const int predicted_class = std::max_element(
        class_to_score.begin(), class_to_score.end(),
        [](auto a, auto b){ return a.second < b.second; }
    )->first;
}

int main() {
    std::string file_name = "data/raw/kaggle_human.txt";
    std::string algorithm = "Smith-Waterman";
    int chuck_len = 64;
    int overlap = 32;
    bool CUDA = false;

    kaggle_pipeline(file_name, algorithm, chuck_len, overlap, CUDA);
}