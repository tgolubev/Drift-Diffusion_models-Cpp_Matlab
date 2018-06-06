#ifndef RECOMBINATION_H
#define RECOMBINATION_H

#include<vector>
#include "parameters.h"
#include "simulation.h"


class Recombo{
public:
    Recombo(Simulation &simul, double k_rec, double E_trap, double n1, double p1)
             : k_rec{k_rec}, E_trap{E_trap}, n1{n1}, p1{p1}, R_Langevin(simul.get_num_cell()){ }       //constructor

    std::vector<double> ComputeR_Langevin(Simulation &simul, const std::vector<double> &n, const std::vector<double> &p);

private:
    double k_rec, E_trap, n1, p1;
    std::vector<double> R_Langevin;
};




#endif // RECOMBINATION_H
