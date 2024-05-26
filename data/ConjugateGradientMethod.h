// ConjugateGradientMethod.h
#ifndef CONJUGATEGRADIENTMETHOD_H
#define CONJUGATEGRADIENTMETHOD_H

#include <vector>
#include "MatrixCalculator.h"

class ConjugateGradientMethod : public MatrixCalculator {
public:
    std::vector<double> conjugate_gradient(const std::vector<std::vector<double>>& Q, const std::vector<double>& b, std::vector<double> x0, double tol=1e-6);
};

#endif // CONJUGATEGRADIENTMETHOD_H