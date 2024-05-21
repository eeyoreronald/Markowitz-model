// ParameterEstimation.cpp
#include "ParameterEstimation.h"
#include <numeric>

ParameterEstimation::ParameterEstimation(const std::vector<std::vector<double>>& returnMatrix) : returnMatrix(returnMatrix) {}

std::vector<double> ParameterEstimation::calculateMeanReturns() {
    std::vector<double> meanReturns;
    for (const auto& assetReturns : returnMatrix) {
        double sum = std::accumulate(assetReturns.begin(), assetReturns.end(), 0.0);
        meanReturns.push_back(sum / assetReturns.size());
    }
    return meanReturns;
}

std::vector<std::vector<double>> ParameterEstimation::calculateCovarianceMatrix() {
    std::vector<double> meanReturns = calculateMeanReturns();
    int numAssets = returnMatrix.size();
    int numReturns = returnMatrix[0].size();

    std::vector<std::vector<double>> covarianceMatrix(numAssets, std::vector<double>(numAssets, 0.0));

    for (int i = 0; i < numAssets; ++i) {
        for (int j = 0; j < numAssets; ++j) {
            double cov = 0.0;
            for (int k = 0; k < numReturns; ++k) {
                cov += (returnMatrix[i][k] - meanReturns[i]) * (returnMatrix[j][k] - meanReturns[j]);
            }
            covarianceMatrix[i][j] = cov / (numReturns - 1);
        }
    }

    return covarianceMatrix;
}