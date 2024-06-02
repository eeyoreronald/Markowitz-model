// In a header file called "MatrixWriter.h"
#ifndef MATRIX_WRITER_H
#define MATRIX_WRITER_H

#include <vector>
#include <string>

void writeMatrixToCSV(const std::vector<std::vector<std::vector<double>>>& matrix, const std::string& filename);

#endif
