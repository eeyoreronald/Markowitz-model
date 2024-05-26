// Quadratic.h
#ifndef QUADRATIC_H
#define QUADRATIC_H

#include <vector>

class Quadratic {
public:
    std::vector<std::vector<double>> createMatrixQ(const std::vector<std::vector<double>>& cov_matrix, const std::vector<double>& mean_returns);
    std::vector<std::vector<double>> calculateBMatrix(const std::vector<double>& target_returns, const std::vector<double>& meanReturns_small);
    std::vector<double> calculateX0(const std::vector<std::vector<double>>& Q, int numberAssets, double initial_guess = 0.5);
    std::vector<std::vector<double>> calculatePortfolioWeights(const std::vector<std::vector<double>>& Q, const std::vector<std::vector<double>>& b_matrix, const std::vector<double>& x0, double tolerance, const std::vector<double>& meanReturns_small);
    std::vector<std::vector<double>> createSmallReturnMatrix(const std::vector<std::vector<double>>& returnMatrix, int numberAssets_small, int numberReturns_small);
    std::vector<std::vector<double>> selectRows(const std::vector<std::vector<double>>& returnMatrix, int start, int end);
    std::vector<std::vector<double>> backtesting(const std::vector<std::vector<double>>& optimal_weights, const std::vector<std::vector<double>>& OOS_returns, const std::vector<double>& target_return);
};

#endif