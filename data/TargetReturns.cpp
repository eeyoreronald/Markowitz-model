// TargetReturns.cpp
#include "TargetReturns.h"

TargetReturns::TargetReturns(double start, double end, int num) {
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        target_returns.push_back(start + i * step);
    }
}

std::vector<double> TargetReturns::getReturns() const {
    return target_returns;
}