#include <iostream>
#include <string>

void pipeline(std::string fileName, ) {
	std::ifstream infile(argv[1]);
	if (!infile) {
		std::cerr << "Error opening file: " << argv[1] << "\n";
		return 1;
	}

	std::string seq1, seq2;

	while (std::getline(infile, seq1) && std::getline(infile, seq2)) {
		SWResult result = smith_waterman(seq1, seq2, 2, -1, -2);

		std::cout << "Score: " << result.score << "\n";
		std::cout << result.aligned_seq1 << "\n"
		          << result.aligned_seq1 << "\n\n";
	}

	return 0;
}
