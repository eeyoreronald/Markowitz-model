//g++ -std=c++17 -c main.cpp
//g++ -std=c++17 -c ParameterEstimation.cpp
//g++ -std=c++17 -c TargetReturns.cpp
//g++ -std=c++17 -c MatrixCalculator.cpp
//g++ -std=c++17 -c ConjugateGradientMethod.cpp
//g++ -std=c++17 -c csv.cpp
//g++ -std=c++17 -c string_to_double.cpp
//g++ -std=c++17 -c readData.cpp
//g++ -std=c++17 -c Quadratic.cpp
//g++ -std=c++17 -c MatrixWriter.cpp


//g++ -std=c++17 -c main.cpp MatrixWriter.cpp Quadratic.cpp readData.cpp string_to_double.cpp ParameterEstimation.cpp TargetReturns.cpp MatrixCalculator.cpp ConjugateGradientMethod.cpp csv.cpp
//g++ -std=c++17 -o portfolioSolver main.o MatrixWriter.o Quadratic.o readData.o csv.o string_to_double.o ParameterEstimation.o TargetReturns.o MatrixCalculator.o ConjugateGradientMethod.o
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
#include "MatrixWriter.h"

using namespace std;

typedef vector<vector<double>> Matrix;  // define a type for a matrix


int main (int argc, char *argv[])
{
    /*
    Matrix A_data = {
    {4, -4},
    {3, 1},
    };
    
    
    // Define the vector b
    vector<double> b_data = { 7,5 };


    vector <double> x0_ed = {0.5, 0.5}; // initial guess for the weights
    cout<<"Start"<<endl;
    ConjugateGradientMethod solver_ed;
    vector<double> results_ed = solver_ed.conjugate_gradient(A_data, b_data, x0_ed);
    cout<<"End"<<endl;
    //print out the result vector
    
    cout<<"\nresults_ed Vector: ";
    for (const auto& x : results_ed) {
        cout << x << " ";
    }
    
    */

    
    int numberAssets=83;    
    int numberReturns=700;
    
    TargetReturns target_1(0.0, 0.1, 21); 
    vector<double> target_returns = target_1.getReturns();  // get the vector of target returns
    double tolerance  = 1e-10;                              // set the tolerance for the conjugate gradient method

    Matrix returnMatrix(numberAssets, vector<double>(numberReturns, 0.0)); // create a matrix to store the returns of the assets
    string fileName="asset_returns.csv";
    readData(returnMatrix,fileName);                        //read the data from the file and store it into the return matrix

    Quadratic portfolio;                                    // create a Quadratic object for portfolio optimization

    vector<vector<vector<double>>> results; 
    vector<vector<vector<double>>> results_2;
    

    for (int i = 0; i < numberReturns - 100; i += 12) {
        int index = i / 12;
        int start = i;
        int mid = i + 100;
        int end = i + 112;
        
        //print out the start, mid and end
        cout<<"\nBacktesting Iteration: "<<index<<endl;
        cout<<"start: "<<start<<" mid: "<<mid<<" end: "<<end<<endl;

        Matrix returnMatrix_IS = portfolio.selectColumns(returnMatrix, start, mid); // select the in-sample data
        Matrix returnMatrix_OOS = portfolio.selectColumns(returnMatrix, mid, end);  // select the out-of-sample data
        
        ParameterEstimation return_IS(returnMatrix_IS);     // create a ParameterEstimation object for the in-sample data

        vector<double> meanReturns_IS = return_IS.calculateMeanReturns(); // calculate the mean returns of the in-sample data
        Matrix covarianceMatrix_IS = return_IS.calculateCovarianceMatrix(); // calculate the covariance matrix of the in-sample data

        //Calculate the matrix Q, b_matrix and x0 for ConjugateGradientMethod
        Matrix Q = portfolio.createMatrixQ(covarianceMatrix_IS, meanReturns_IS); // create the Q matrix
        Matrix b_matrix = portfolio.calculateBMatrix(target_returns, meanReturns_IS); // create the b matrix
        vector <double> x0 = portfolio.calculateX0(Q, numberAssets, 0.5);       // calculate the initial guess for the weights




        Matrix df_weights = portfolio.calculatePortfolioWeights(Q, b_matrix, x0, tolerance, meanReturns_IS); // calculate the optimal weights
        Matrix OOS_returns = portfolio.backtesting(df_weights, returnMatrix_OOS, target_returns); // backtesting the optimal weights

        //store the results
        results.push_back(OOS_returns);
        results_2.push_back(df_weights);
    }
    
    //write the results to a csv file
    writeMatrixToCSV(results, "OOS_returns.csv");
    writeMatrixToCSV(results_2, "df_weights.csv");

    
    
    return 0;

}



