// Quadratic.cpp
#include "Quadratic.h"
#include "ConjugateGradientMethod.h"
#include "ParameterEstimation.h"
#include "MatrixCalculator.h"

#include <iostream>

std::vector<std::vector<double>> Quadratic::createMatrixQ(const std::vector<std::vector<double>>& cov_matrix, const std::vector<double>& mean_returns) {
    int len = mean_returns.size();
    std::vector<std::vector<double>> Q(len + 2, std::vector<double>(len + 2, 0));

    // Copy cov_matrix into Q
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            Q[i][j] = cov_matrix[i][j];
        }
    }

    // Copy mean_returns and -1 into Q
    for (int i = 0; i < len; i++) {
        Q[len][i] = -mean_returns[i];
        Q[len + 1][i] = -1;
    }

    // Copy array_ and array_2 into Q
    for (int i = 0; i < len; i++) {
        Q[i][len] = -mean_returns[i];
        Q[i][len + 1] = -1;
    }

    // Print Q
    std::cout << "\n_______\n Matrix Q:\n";
    for (const auto& row : Q) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << "\n";
    }
    std::cout << "(" << Q.size() << ", " << Q[0].size() << ")\n";

    return Q;
}

std::vector<std::vector<double>> Quadratic::calculateBMatrix(const std::vector<double>& target_returns, const std::vector<double>& meanReturns_small) {
    int n = target_returns.size();
    int mean_returns_len = meanReturns_small.size();

    std::vector<std::vector<double>> b_matrix(n, std::vector<double>(mean_returns_len + 2, 0.0));

    for (int i = 0; i < n; ++i) {
        b_matrix[i][mean_returns_len] = -target_returns[i];
        b_matrix[i][mean_returns_len + 1] = -1.0;
    }

    //print out b_matrix
    std::cout << "\n_______\nPrinting out b_matrix\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < mean_returns_len + 2; j++) {
            std::cout << b_matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
    return b_matrix;
}

std::vector<double> Quadratic::calculateX0(const std::vector<std::vector<double>>& Q, int numberAssets, double initial_guess) {
    int x0_size = Q.size();
    double fillValue = initial_guess;
    std::vector<double> x0(x0_size, fillValue);

    //print out x0
    std::cout << "\n_______\nPrinting out x0 in the main function\n";
    for(int i = 0; i < numberAssets + 2; i++) {
        std::cout << x0[i] << " ";
    }

    return x0;
}

std::vector<std::vector<double>> Quadratic::calculatePortfolioWeights(const std::vector<std::vector<double>>& Q, const std::vector<std::vector<double>>& b_matrix, const std::vector<double>& x0, double tolerance, const std::vector<double>& meanReturns_small) {
    std::vector<std::vector<double>> weights_list;

    ConjugateGradientMethod solver;

    for (int i = 0; i < b_matrix.size(); ++i) {
        std::vector<double> weight = solver.conjugate_gradient(Q, b_matrix[i], x0, tolerance);
        weight.resize(meanReturns_small.size());
        weights_list.push_back(weight);
    }

    std::cout << "\nWeights Matrix: \n";
    for (const auto& weights : weights_list) {
        for (const auto& x : weights) {
            std::cout << x << " ";
        }
        std::cout << "\n";
    }

    return weights_list;
}

std::vector<std::vector<double>> Quadratic::createSmallReturnMatrix(const std::vector<std::vector<double>>& returnMatrix, int numberAssets_small, int numberReturns_small) {
    std::vector<std::vector<double>> returnMatrix_small(numberAssets_small, std::vector<double>(numberReturns_small, 0.0));

    for(int i = 0; i < numberAssets_small; i++) {
        for(int j = 0; j < numberReturns_small; j++) {
            returnMatrix_small[i][j] = returnMatrix[i][j];
        }
    }

    return returnMatrix_small;
}

std::vector<std::vector<double>> Quadratic::selectRows(const std::vector<std::vector<double>>& returnMatrix, int start, int end) {
    std::vector<std::vector<double>> subset;
    for (int i = start; i < end; ++i) {
        subset.push_back(returnMatrix[i]);
    }
    return subset;
}

std::vector<std::vector<double>> Quadratic::backtesting(const std::vector<std::vector<double>>& optimal_weights, const std::vector<std::vector<double>>& OOS_returns, const std::vector<double>& target_return) {
    ParameterEstimation return_OOS(OOS_returns);
    std::vector<double> mean_returns_OOS = return_OOS.calculateMeanReturns();
    std::vector<std::vector<double>> covariance_martix_OOS = return_OOS.calculateCovarianceMatrix();

    MatrixCalculator calc;
    std::vector<std::vector<double>> res;

    for (int index = 0; index < optimal_weights.size(); ++index) {
        std::vector<double> weights = optimal_weights[index];
        double targ_ret = target_return[index];

        double act_ave_return = calc.dotProduct(mean_returns_OOS, weights);

        std::vector<double> temp = calc.matrixVectorProduct(covariance_martix_OOS, weights);
        double pf_cov = calc.dotProduct(weights, temp);

        res.push_back({targ_ret, act_ave_return, pf_cov});
    }

    return res;
}