// MatrixCalculator.cpp
#include "MatrixCalculator.h"

double MatrixCalculator::dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0;
    for (size_t i = 0; i < a.size(); i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

std::vector<double> MatrixCalculator::matrixVectorProduct(const std::vector<std::vector<double>>& a, const std::vector<double>& x) {
    std::vector<double> result(a.size(), 0);
    for (size_t i = 0; i < a.size(); i++) {
        for (size_t j = 0; j < a[0].size(); j++) {
            result[i] += a[i][j] * x[j];
        }
    }
    return result;
}

std::vector<double> MatrixCalculator::vectorSubtraction(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<double> MatrixCalculator::vectorAddition(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<double> MatrixCalculator::scalarVectorProduct(double scalar, const std::vector<double>& v) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        result[i] = scalar * v[i];
    }
    return result;
}