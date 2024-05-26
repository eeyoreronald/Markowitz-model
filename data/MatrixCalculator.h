// MatrixCalculator.h
#ifndef MATRIXCALCULATOR_H
#define MATRIXCALCULATOR_H

#include <vector>

class MatrixCalculator {
public:
    double dotProduct(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> matrixVectorProduct(const std::vector<std::vector<double>>& a, const std::vector<double>& x);
    std::vector<double> vectorSubtraction(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> vectorAddition(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> scalarVectorProduct(double scalar, const std::vector<double>& v);
};

#endif // MATRIXCALCULATOR_H