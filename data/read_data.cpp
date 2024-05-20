//g++ -std=c++17 -c read_data.cpp
//g++ -std=c++17 -c csv.cpp
//g++ -std=c++17 -o portfolioSolver csv.o read_data.o
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

using namespace std;


double string_to_double( const string& s );
void readData(vector<vector<double> >& data, const string& fileName);


//---------------------- Class implmentation -----------------------------------

class ParameterEstimation {
private:
    vector<vector<double> > returnMatrix;

public:
    ParameterEstimation(const vector<vector<double> >& returnMatrix) : returnMatrix(returnMatrix) {}

    vector<double> calculateMeanReturns() {
        vector<double> meanReturns;
        for (const auto& assetReturns : returnMatrix) {
            double sum = accumulate(assetReturns.begin(), assetReturns.end(), 0.0);
            meanReturns.push_back(sum / assetReturns.size());
        }
        return meanReturns;
    }

    vector<vector<double> > calculateCovarianceMatrix() {
        vector<double> meanReturns = calculateMeanReturns();
        int numAssets = returnMatrix.size();
        int numReturns = returnMatrix[0].size();

        vector<vector<double> > covarianceMatrix(numAssets, vector<double>(numAssets, 0.0));

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


};

class TargetReturns {
private:
    vector<double> target_returns;

public:
    TargetReturns(double start, double end, int num) {
        double step = (end - start) / (num - 1);
        for (int i = 0; i < num; ++i) {
            target_returns.push_back(start + i * step);
        }
    }

    vector<double> getReturns() const {
        return target_returns;
    }
};

class MatrixCalculator {
public:
    double dotProduct(const vector<double>& a, const vector<double>& b) {
        double sum = 0;
        for (size_t i = 0; i < a.size(); i++) {
            sum += a[i] * b[i];
        }
        return sum;
    }

    vector<double> matrixVectorProduct(const vector<vector<double> >& a, const vector<double>& x) {
        vector<double> result(a.size(), 0);
        for (size_t i = 0; i < a.size(); i++) {
            for (size_t j = 0; j < a[0].size(); j++) {
                result[i] += a[i][j] * x[j];
            }
        }
        return result;
    }

    vector<double> vectorSubtraction(const vector<double>& a, const vector<double>& b) {
        vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); i++) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    vector<double> vectorAddition(const vector<double>& a, const vector<double>& b) {
        vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); i++) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    vector<double> scalarVectorProduct(double scalar, const vector<double>& v) {
        vector<double> result(v.size());
        for (size_t i = 0; i < v.size(); i++) {
            result[i] = scalar * v[i];
        }
        return result;
    }
};

class ConjugateGradientMethod : public MatrixCalculator {
public:
    vector<double> conjugate_gradient(const vector<vector<double> >& Q, const vector<double>& b, vector<double> x0, double tol=1e-6) {
        vector<double> s0 = vectorSubtraction(b, matrixVectorProduct(Q, x0));
        vector<double> p0 = s0;
        int i = 0;
        while (dotProduct(s0, s0) > tol) {
            double alpha = dotProduct(s0, s0) / dotProduct(matrixVectorProduct(Q, p0), p0);
            x0 = vectorAddition(x0, scalarVectorProduct(alpha, p0));
            vector<double> s1 = vectorSubtraction(s0, scalarVectorProduct(alpha, matrixVectorProduct(Q, p0)));
            double beta = dotProduct(s1, s1) / dotProduct(s0, s0);
            p0 = vectorAddition(s1, scalarVectorProduct(beta, p0));
            s0 = s1;
            i++;
            cout << "This is iteration number: " << i << "\n";
            cout << "The current x0 is: ";
            for (const auto& x : x0) {
                cout << x << " ";
            }
            cout << "\n";
        }
        return x0;
    }
};

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
            Q[len][i] = mean_returns[i];
            Q[len + 1][i] = -1;
        }

        // Copy array_ and array_2 into Q
        for (int i = 0; i < len; i++) {
            Q[i][len] = mean_returns[i];
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

    Quadratic quadratic;
    vector<vector<double>> returnMatrix_small = quadratic.createSmallReturnMatrix(returnMatrix, numberAssets_small, numberReturns_small);
    

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
    vector<vector<double> > Q = quadratic.createMatrixQ(covarianceMatrix_small, meanReturns_small);
    vector<vector<double> > b_matrix = quadratic.calculateBMatrix(target_returns, meanReturns_small);
    vector <double> x0 = quadratic.calculateX0(Q, numberAssets_small, 0.5);
    
    //For large dataset
    //Q = quadratic.createMatrixQ(covarianceMatrix_large, meanReturns_large);
    //b_matrix = quadratic.calculateBMatrix(target_returns, meanReturns_large);
    //x0 = quadratic.calculateX0(Q, numberAssets, 0.5);


    //ConjugateGradientMethod solver, for portfolio optimization
    double tolerance  = 1e-10;
    vector<vector<double> > weights_matrix_small = quadratic.calculatePortfolioWeights(Q, b_matrix, x0, tolerance, meanReturns_small);

    
    //Backtesting in the following code



    return 0;
}




//--------------------------------------------------------------------------------
// Original code for readData function
//--------------------------------------------------------------------------------



double string_to_double( const string& s )
{
    istringstream i(s);
    double x;
    if (!(i >> x))
        return 0;
    return x;
} 

void readData(vector<vector<double> >& data, const string& fileName)
{
    ifstream file(fileName);
    Csv csv(file);
    string line;
    if (file.is_open())
    {
        int i = 0;
        while (csv.getline(line) != 0) {
            cout << "Enter while loop, i=" << i << "\n";
            for (int j = 0; j < csv.getnfield(); j++)
            {
                double temp = string_to_double(csv.getfield(j));
                cout << "Asset " << j << ", Return " << i << "=" << temp << "\n";
                cout << "j: " << j << ", i: " << i << "\n";
                data[j][i] = temp;
            }
            i++;
            cout << "end while loop, i=" << i << "\n";
        }
        cout << "Read successfully" << endl;
        file.close();
        cout << "File closed" << endl;


    }
    else {
        cout << fileName << " missing\n";
        exit(0);
    }
    return;
}

