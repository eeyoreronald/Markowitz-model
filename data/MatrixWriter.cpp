// In a cpp file called "MatrixWriter.cpp"
#include "MatrixWriter.h"
#include <iostream>
#include <fstream>

void writeMatrixToCSV(const std::vector<std::vector<std::vector<double>>>& matrix, const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cout << "Failed to open the file." << std::endl;
        return;
    }

    for (const auto& subMatrix : matrix) {
        for (const auto& row : subMatrix) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << row[i];
                if (i != row.size() - 1) {
                    file << ",";
                }
            }
            file << "\n";
        }
        file << "\n";
    }

    file.close();
}