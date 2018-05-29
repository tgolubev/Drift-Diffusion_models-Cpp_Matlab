#include<vector>
#include "parameters.h"
#include<iostream>

void ComputeR_Langevin(const std::vector<double> &n, const std::vector<double> &p, std::vector<double> &R_Langevin){
    for(int i = 1;i<= num_cell-1;i++){
        R_Langevin[i] = k_rec*(Nsqrd*n[i]*p[i] - n1*p1);
        if(R_Langevin[i] < 0.0)  R_Langevin[i] = 0.0;  //negative recombo is unphysical
    }
}

