// ConjugateGradientMethod.cpp
#include "ConjugateGradientMethod.h"

std::vector<double> ConjugateGradientMethod::conjugate_gradient(const std::vector<std::vector<double>>& Q, const std::vector<double>& b, std::vector<double> x0, double tol) {
    std::vector<double> s0 = vectorSubtraction(b, matrixVectorProduct(Q, x0));
    std::vector<double> p0 = s0;
    int i = 0;
    while (dotProduct(s0, s0) > tol) {
        double alpha = dotProduct(s0, s0) / dotProduct(matrixVectorProduct(Q, p0), p0);
        x0 = vectorAddition(x0, scalarVectorProduct(alpha, p0));
        std::vector<double> s1 = vectorSubtraction(s0, scalarVectorProduct(alpha, matrixVectorProduct(Q, p0)));
        double beta = dotProduct(s1, s1) / dotProduct(s0, s0);
        p0 = vectorAddition(s1, scalarVectorProduct(beta, p0));
        s0 = s1;
        i++;
    }
    return x0;
}