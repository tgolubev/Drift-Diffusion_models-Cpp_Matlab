#ifndef RECOMBINATION_H
#define RECOMBINATION_H

#include<vector>
#include "constants.h"
#include "parameters.h"


class Recombo{
public:
    Recombo(Parameters &params, double k_rec)
             : k_rec{k_rec}, R_Langevin(params.num_cell)
    {
        E_trap = params.active_VB + params.E_gap/2.0;  //trap assisted recombo is most effective when trap is located mid-gap--> take VB and add 1/2 of bandgap
        n1 = params.N_LUMO*exp(-(params.active_CB - E_trap)/Vt);
        p1 = params.N_HOMO*exp(-(E_trap - params.active_VB)/Vt);       //constructor
    }

    std::vector<double> ComputeR_Langevin(Parameters &params, const std::vector<double> &n, const std::vector<double> &p);
    void set_k_rec(double k_rec_input){ k_rec = k_rec_input;}
private:
    double k_rec;
    double E_trap;
    double n1;
    double p1;
    std::vector<double> R_Langevin;
};


#endif // RECOMBINATION_H
