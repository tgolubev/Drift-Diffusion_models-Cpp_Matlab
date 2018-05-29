#ifndef RECOMBINATION_H
#define RECOMBINATION_H

#include<vector>
#include "parameters.h"

void ComputeR_Langevin(const std::vector<double> &n, const std::vector<double> &p, std::vector<double> &R_Langevin);


#endif // RECOMBINATION_H
