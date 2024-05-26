//g++ -std=c++17 -c main.cpp
//g++ -std=c++17 -c ParameterEstimation.cpp
//g++ -std=c++17 -c TargetReturns.cpp
//g++ -std=c++17 -c MatrixCalculator.cpp
//g++ -std=c++17 -c ConjugateGradientMethod.cpp
//g++ -std=c++17 -c csv.cpp
//g++ -std=c++17 -c string_to_double.cpp
//g++ -std=c++17 -c readData.cpp
//g++ -std=c++17 -c Quadratic.cpp


//g++ -std=c++17 -c main.cpp Quadratic.cpp readData.cpp string_to_double.cpp ParameterEstimation.cpp TargetReturns.cpp MatrixCalculator.cpp ConjugateGradientMethod.cpp csv.cpp
//g++ -std=c++17 -o portfolioSolver main.o Quadratic.o readData.o csv.o string_to_double.o ParameterEstimation.o TargetReturns.o MatrixCalculator.o ConjugateGradientMethod.o
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
#include "Quadratic.h"

using namespace std;



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
