#ifndef BERNOULLI_H
#define BERNOULLI_H


#include<vector>
#include "parameters.h"

void BernoulliFnc_n(const std::vector<double> &V, std::vector<double> &B_n1, std::vector<double> &B_n2);
void BernoulliFnc_p(const std::vector<double> &V, std::vector<double> &B_p1, std::vector<double> &B_p2);

#endif // BERNOULLI_H
