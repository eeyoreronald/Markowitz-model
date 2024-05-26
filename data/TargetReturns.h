// TargetReturns.h
#ifndef TARGETRETURNS_H
#define TARGETRETURNS_H

#include <vector>

class TargetReturns {
private:
    std::vector<double> target_returns;

public:
    TargetReturns(double start, double end, int num);
    std::vector<double> getReturns() const;
};

#endif // TARGETRETURNS_H