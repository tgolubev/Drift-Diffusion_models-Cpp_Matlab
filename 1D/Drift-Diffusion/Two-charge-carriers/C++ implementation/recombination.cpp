#include<vector>
#include<iostream>
#include "recombination.h"

std::vector<double> Recombo::ComputeR_Langevin(Simulation &simul, const std::vector<double> &n, const std::vector<double> &p){
    for(int i = 1;i<= simul.get_num_cell()-1;i++){
        R_Langevin[i] = k_rec*(N*N*n[i]*p[i] - n1*p1);
        if(R_Langevin[i] < 0.0)  R_Langevin[i] = 0.0;  //negative recombo is unphysical
    }

    return R_Langevin;
}

