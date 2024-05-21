// ParameterEstimation.h
#ifndef PARAMETERESTIMATION_H
#define PARAMETERESTIMATION_H

#include <vector>

class ParameterEstimation {
private:
    std::vector<std::vector<double>> returnMatrix;

public:
    ParameterEstimation(const std::vector<std::vector<double>>& returnMatrix);

    std::vector<double> calculateMeanReturns();
    std::vector<std::vector<double>> calculateCovarianceMatrix();
};

#endif // PARAMETERESTIMATION_H