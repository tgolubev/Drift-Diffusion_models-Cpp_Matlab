#ifndef RECOMBINATION_H
#define RECOMBINATION_H

#include<vector>
#include "parameters.h"
#include "simulation.h"


class Recombo{
public:
    Recombo(Simulation &simul, double k_rec)
             : k_rec{k_rec}, E_trap{E_trap}, n1{n1}, p1{p1}, R_Langevin(simul.get_num_cell()){ }       //constructor

    std::vector<double> ComputeR_Langevin(const std::vector<double> &n, const std::vector<double> &p);

private:
    double k_rec;
    const double E_trap = active_VB + E_gap/2.0;  //this is assumption used--> trap assisted recombo is most effective when trap is located mid-gap--> take VB and add 1/2 of bandgap
    const double n1 = N_LUMO*exp(-(active_CB - E_trap)/Vt);
    const double p1 = N_HOMO*exp(-(E_trap - active_VB)/Vt);
    std::vector<double> R_Langevin;
};




#endif // RECOMBINATION_H
