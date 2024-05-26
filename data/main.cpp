//g++ -std=c++17 -c main.cpp
//g++ -std=c++17 -c ParameterEstimation.cpp
//g++ -std=c++17 -c TargetReturns.cpp
//g++ -std=c++17 -c MatrixCalculator.cpp
//g++ -std=c++17 -c ConjugateGradientMethod.cpp
//g++ -std=c++17 -c csv.cpp
//g++ -std=c++17 -c string_to_double.cpp
//g++ -std=c++17 -c readData.cpp


//g++ -std=c++17 -c main.cpp readData.cpp string_to_double.cpp ParameterEstimation.cpp TargetReturns.cpp MatrixCalculator.cpp ConjugateGradientMethod.cpp csv.cpp
//g++ -std=c++17 -o portfolioSolver main.o readData.o csv.o string_to_double.o ParameterEstimation.o TargetReturns.o MatrixCalculator.o ConjugateGradientMethod.o
//./portfolioSolver

#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>

#include "csv.h"
#include "ParameterEstimation.h"
#include "TargetReturns.h"
#include "MatrixCalculator.h"
#include "ConjugateGradientMethod.h"
#include "string_to_double.h"
#include "readData.h"

using namespace std;


//---------------------- Class implmentation -----------------------------------

class Quadratic {
public:
    vector<vector<double> > createMatrixQ(const vector<vector<double> >& cov_matrix, const vector<double>& mean_returns) {
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
    vector<vector<double> > calculateBMatrix(const vector<double>& target_returns, const vector<double>& meanReturns_small) {
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
    vector<double> calculateX0(const vector<vector<double>>& Q, int numberAssets, double initial_guess = 0.5) {
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
    vector<vector<double> > calculatePortfolioWeights(const vector<vector<double>>& Q, const vector<vector<double>>& b_matrix, const vector<double>& x0, double tolerance, const vector<double>& meanReturns_small) {
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

    vector<vector<double> > createSmallReturnMatrix(const vector<vector<double>>& returnMatrix, int numberAssets_small, int numberReturns_small) {
        vector<vector<double> > returnMatrix_small(numberAssets_small, vector<double>(numberReturns_small, 0.0));

        for(int i = 0; i < numberAssets_small; i++) {
            for(int j = 0; j < numberReturns_small; j++) {
                returnMatrix_small[i][j] = returnMatrix[i][j];
            }
        }

        return returnMatrix_small;
    }

    vector<vector<double>> selectRows(const vector<vector<double>>& returnMatrix, int start, int end) {
        vector<vector<double>> subset;
        for (int i = start; i < end; ++i) {
            subset.push_back(returnMatrix[i]);
        }
        return subset;
    }

    vector<vector<double>> backtesting (vector<vector<double>> optimal_weights, vector<vector<double>> OOS_returns, vector<double> target_return){
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

};

//---------------------- Class implmentation -----------------------------------



int main (int argc, char *argv[])
{
    
    int numberAssets=83;    
    int numberReturns=700;

    vector<vector<double> > returnMatrix(numberAssets, vector<double>(numberReturns, 0.0));
    
    string fileName="asset_returns.csv";
    readData(returnMatrix,fileName);
    

    //read the data from the file and store it into the return matrix



    //Create a matrix to store the first 5 assets and the first 10 returns
    int numberAssets_small=5;
    int numberReturns_small=10;
    Quadratic portfolio;
    vector<vector<double>> returnMatrix_small = portfolio.createSmallReturnMatrix(returnMatrix, numberAssets_small, numberReturns_small);
    

    //--------------------------------------------------------------------------------
    //Remark to myself
    //i = asset number
    //j = return number
    //returnMatrix[i][j] stores the asset i, return j value
    //--------------------------------------------------------------------------------

    //Create a ParameterEstimation object
    ParameterEstimation return_small(returnMatrix_small);
    ParameterEstimation return_large(returnMatrix);


    vector<double> meanReturns_small = return_small.calculateMeanReturns();
    vector<vector<double> > covarianceMatrix_small = return_small.calculateCovarianceMatrix();
    
    vector<double> meanReturns_large = return_large.calculateMeanReturns();
    vector<vector<double> > covarianceMatrix_large = return_large.calculateCovarianceMatrix();


   //Vector to store the target returns
    TargetReturns target_1(0.0, 0.1, 20);
    vector<double> target_returns = target_1.getReturns();

    //Calculate the matrix Q, b_matrix and x0 for ConjugateGradientMethod
    vector<vector<double> > Q = portfolio.createMatrixQ(covarianceMatrix_small, meanReturns_small);
    vector<vector<double> > b_matrix = portfolio.calculateBMatrix(target_returns, meanReturns_small);
    vector <double> x0 = portfolio.calculateX0(Q, numberAssets_small, 0.5);
    
    //For large dataset
    //Q = portfolio.createMatrixQ(covarianceMatrix_large, meanReturns_large);
    //b_matrix = portfolio.calculateBMatrix(target_returns, meanReturns_large);
    //x0 = portfolio.calculateX0(Q, numberAssets, 0.5);


    //ConjugateGradientMethod solver, for portfolio optimization
    double tolerance  = 1e-10;
    vector<vector<double> > weights_matrix_small = portfolio.calculatePortfolioWeights(Q, b_matrix, x0, tolerance, meanReturns_small);
    //vector<vector<double> > weights_matrix_large = portfolio.calculatePortfolioWeights(Q, b_matrix, x0, tolerance, meanReturns_large);

    
    //Backtesting in the following code
    //vector<vector<double>> backtesting (vector<vector<double>> optimal_weights, vector<vector<double>> OOS_returns, vector<double> target_return)

    vector<vector<double>> OOS_return = portfolio.backtesting(weights_matrix_small, returnMatrix_small, target_returns);
    
    //print out OOS_return
    cout << "\n_______\nPrinting out OOS_return\n";
    for(int i = 0; i < OOS_return.size(); i++) {
        for(int j = 0; j < OOS_return[i].size(); j++) {
            cout << OOS_return[i][j] << " ";
        }
        cout << "\n";
    }

    //note to myself
    //OOS_return[i][0] = target return
    //OOS_return[i][1] = actual average return
    //OOS_return[i][2] = portfolio covariance

    //-------------------------------------------------------------------------------- Backtesting implmentations
    
    /*
    //vector<vector<double> > returnMatrix(numberAssets, vector<double>(numberReturns, 0.0));
    //string fileName="asset_returns.csv";
    //readData(returnMatrix,fileName);
    //int numberReturns = 700;
    int N = (numberReturns - 100) / 12;

    
    TargetReturns target_1(0.0, 0.1, 20);
    vector<double> target_returns = target_1.getReturns();
    double tolerance  = 1e-10;
    
    if (numberReturns > returnMatrix.size()) {
        throw length_error("numberReturns is greater than the number of rows in returnMatrix");
    }
    else
        cout<<"numberReturns is less than the number of rows in returnMatrix"<<endl;



    for (int i = 0; i < numberReturns - 100; i += 12) {
        int index = i / 12;
        int start = i;
        int mid = i + 100;
        int end = i + 112;
        vector<vector<double>> returnMatrix_IS = portfolio.selectRows(returnMatrix, start, mid);
        vector<vector<double>> returnMatrix_OOS = portfolio.selectRows(returnMatrix, mid, end);
        ParameterEstimation return_IS(returnMatrix_IS);


        
        vector<double> meanReturns_IS = return_IS.calculateMeanReturns();
        vector<vector<double> > covarianceMatrix_IS = return_IS.calculateCovarianceMatrix();
        
        Q = portfolio.createMatrixQ(covarianceMatrix_IS, meanReturns_IS);
        b_matrix = portfolio.calculateBMatrix(target_returns, meanReturns_IS);
        x0 = portfolio.calculateX0(Q, numberAssets, 0.5);

        vector<vector<double> > weights_matrix_IS = portfolio.calculatePortfolioWeights(Q, b_matrix, x0, tolerance, meanReturns_IS);
        vector<vector<double>> OOS_return = portfolio.backtesting(weights_matrix_IS, returnMatrix_OOS, target_returns);
        cout << "index: " << index << "\n";
        


    }
    
    */
    
    
    
    return 0;
    

}
