#include "Quadratic.h"
#include <iostream>

// Implementations of the Quadratic class methods go here
// For example:

vector<vector<double> > Quadratic::createMatrixQ(const vector<vector<double> >& cov_matrix, const vector<double>& mean_returns) {
    // method implementation goes here
}

vector<vector<double> > Quadratic::createMatrixQ(const vector<vector<double> >& cov_matrix, const vector<double>& mean_returns) {
        int len = mean_returns.size();
        vector<vector<double> > Q(len + 2, vector<double>(len + 2, 0));

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
        cout << "\n_______\n Matrix Q:\n";
        for (const auto& row : Q) {
            for (const auto& elem : row) {
                cout << elem << " ";
            }
            cout << "\n";
        }
        cout << "(" << Q.size() << ", " << Q[0].size() << ")\n";

        return Q;
    }

vector<vector<double> > Quadratic::calculateBMatrix(const vector<double>& target_returns, const vector<double>& meanReturns_small) {
    int n = target_returns.size();
    int mean_returns_len = meanReturns_small.size();

    vector<vector<double>> b_matrix(n, vector<double>(mean_returns_len + 2, 0.0));

    for (int i = 0; i < n; ++i) {
        b_matrix[i][mean_returns_len] = -target_returns[i];
        b_matrix[i][mean_returns_len + 1] = -1.0;
    }

    //print out b_matrix
    cout << "\n_______\nPrinting out b_matrix\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < mean_returns_len + 2; j++) {
            cout << b_matrix[i][j] << " ";
        }
        cout << "\n";
    }
    return b_matrix;
}
vector<double> Quadratic::calculateX0(const vector<vector<double>>& Q, int numberAssets, double initial_guess = 0.5) {
    int x0_size = Q.size();
    double fillValue = initial_guess;
    vector<double> x0(x0_size, fillValue);

    //print out x0
    cout << "\n_______\nPrinting out x0 in the main function\n";
    for(int i = 0; i < numberAssets + 2; i++) {
        cout << x0[i] << " ";
    }

    return x0;
}
vector<vector<double> > Quadratic::calculatePortfolioWeights(const vector<vector<double>>& Q, const vector<vector<double>>& b_matrix, const vector<double>& x0, double tolerance, const vector<double>& meanReturns_small) {
    vector<vector<double> > weights_list;

    ConjugateGradientMethod solver;

    for (int i = 0; i < b_matrix.size(); ++i) {
        vector<double> weight = solver.conjugate_gradient(Q, b_matrix[i], x0, tolerance);
        weight.resize(meanReturns_small.size());
        weights_list.push_back(weight);
    }

    cout << "\nWeights Matrix: \n";
    for (const auto& weights : weights_list) {
        for (const auto& x : weights) {
            cout << x << " ";
        }
        cout << "\n";
    }

    return weights_list;
}

vector<vector<double> > Quadratic::createSmallReturnMatrix(const vector<vector<double>>& returnMatrix, int numberAssets_small, int numberReturns_small) {
    vector<vector<double> > returnMatrix_small(numberAssets_small, vector<double>(numberReturns_small, 0.0));

    for(int i = 0; i < numberAssets_small; i++) {
        for(int j = 0; j < numberReturns_small; j++) {
            returnMatrix_small[i][j] = returnMatrix[i][j];
        }
    }

    return returnMatrix_small;
}

vector<vector<double>> Quadratic::selectRows(const vector<vector<double>>& returnMatrix, int start, int end) {
    vector<vector<double>> subset;
    for (int i = start; i < end; ++i) {
        subset.push_back(returnMatrix[i]);
    }
    return subset;
}

vector<vector<double>> Quadratic::backtesting (vector<vector<double>> optimal_weights, vector<vector<double>> OOS_returns, vector<double> target_return){
    ParameterEstimation return_OOS(OOS_returns);
    vector <vector<double>> covariance_martix_OOS = return_OOS.calculateCovarianceMatrix();
    vector <double> mean_returns_OOS = return_OOS.calculateMeanReturns();

    MatrixCalculator calc;
    vector<vector<double>> res;

    for (int index = 0; index < optimal_weights.size(); ++index) {
        vector<double> weights = optimal_weights[index];
        double targ_ret = target_return[index];

        double act_ave_return = calc.dotProduct(mean_returns_OOS, weights);

        vector<double> temp = calc.matrixVectorProduct(covariance_martix_OOS, weights);
        double pf_cov = calc.dotProduct(weights, temp);

        res.push_back({targ_ret, act_ave_return, pf_cov});
    }

    return res;

}