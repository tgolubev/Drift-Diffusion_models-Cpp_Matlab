#ifndef BERNOULLI_H
#define BERNOULLI_H


#include<vector>
#include "parameters.h"
#include "simulation.h"  //needs this to know what Simulation is

void BernoulliFnc_n(Simulation &simul, const std::vector<double> &V, std::vector<double> &B_n1, std::vector<double> &B_n2);
void BernoulliFnc_p(Simulation &simul, const std::vector<double> &V, std::vector<double> &B_p1, std::vector<double> &B_p2);

#endif // BERNOULLI_H
