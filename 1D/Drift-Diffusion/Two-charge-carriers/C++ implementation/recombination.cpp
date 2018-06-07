#include<vector>
#include<iostream>
#include "recombination.h"

std::vector<double> Recombo::ComputeR_Langevin(Parameters &params, const std::vector<double> &n, const std::vector<double> &p){
    for(int i = 1;i<p.size();i++){
        R_Langevin[i] = k_rec*(params.N*params.N*n[i]*p[i] - n1*p1);
        if(R_Langevin[i] < 0.0)  R_Langevin[i] = 0.0;  //negative recombo is unphysical
    }

    return R_Langevin;
}

