#include <vector>
using namespace std;

class Quadratic {
public:
    vector<vector<double> > createMatrixQ(const vector<vector<double> >& cov_matrix, const vector<double>& mean_returns);
    vector<vector<double> > calculateBMatrix(const vector<double>& target_returns, const vector<double>& meanReturns_small);
    vector<double> calculateX0(const vector<vector<double>>& Q, int numberAssets, double initial_guess = 0.5);
    vector<vector<double> > calculatePortfolioWeights(const vector<vector<double>>& Q, const vector<vector<double>>& b_matrix, const vector<double>& x0, double tolerance, const vector<double>& meanReturns_small);
    vector<vector<double> > createSmallReturnMatrix(const vector<vector<double>>& returnMatrix, int numberAssets_small, int numberReturns_small);
    vector<vector<double>> selectRows(const vector<vector<double>>& returnMatrix, int start, int end);
    vector<vector<double>> backtesting (vector<vector<double>> optimal_weights, vector<vector<double>> OOS_returns, vector<double> target_return);
};