#include<vector>
#include<iostream>
#include "recombination.h"

Recombo:: Recombo(const Parameters &params)      //constructor
{
    k_rec = params.k_rec;
    R_Langevin.resize(params.num_cell);
    E_trap = params.active_VB + params.E_gap/2.0;  //trap assisted recombo is most effective when trap is located mid-gap--> take VB and add 1/2 of bandgap
    n1 = params.N_LUMO*exp(-(params.active_CB - E_trap)/Vt);
    p1 = params.N_HOMO*exp(-(E_trap - params.active_VB)/Vt);
}

std::vector<double> Recombo::ComputeR_Langevin(const Parameters &params, const std::vector<double> &n, const std::vector<double> &p)
{
    for (int i = 1; i < p.size(); i++) {
        R_Langevin[i] = k_rec*(params.N_dos*params.N_dos*n[i]*p[i] - n1*p1);
        if (R_Langevin[i] < 0.0)  R_Langevin[i] = 0.0;  //negative recombo is unphysical
    }

    //for now just use n and p vectors....

    return R_Langevin;
}

