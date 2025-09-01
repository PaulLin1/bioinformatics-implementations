#include <iostream>
#include <fstream>
#include <string>

void analyze_fasta_file(const std::string& file_path) {
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << file_path << std::endl;
        return;
    }

    std::string line;
    std::string seq;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!seq.empty()) {
                std::cout << seq << std::endl;
                seq.clear();
            }
        } else {
            seq += line;
        }
    }
}

int main() {
    const std::string FASTA_FILE_PATH = "/mnt/scratch/linpaul1/1uniprot_trembl.fasta";
    analyze_fasta_file(FASTA_FILE_PATH);
    return 0;
}
