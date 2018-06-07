#include <vector>
#include <cmath>
#include <iostream>
#include "bernoulli.h"


void BernoulliFnc_n(Simulation &simul, const std::vector<double> &V, std::vector<double> &B_n1, std::vector<double> &B_n2){  //NOTE: this V has the right side bndry condition appended
    std::vector<double> dV(simul.get_num_cell()+1);

    for(int i = 1; i<=simul.get_num_cell();i++){
        dV[i] =  V[i]-V[i-1];
    }

    for(int i = 1; i<=simul.get_num_cell();i++){
        B_n1[i] = dV[i]/(exp(dV[i])-1.0);
        B_n2[i] = B_n1[i]*exp(dV[i]);
    }
}


void BernoulliFnc_p(Simulation &simul, const std::vector<double> &V, std::vector<double> &B_p1, std::vector<double> &B_p2){  //NOTE: this V has the right side bndry condition appended
    std::vector<double> dV(simul.get_num_cell()+1);

    for(int i = 1; i<=simul.get_num_cell();i++){
        dV[i] =  V[i]-V[i-1];
    }

    for(int i = 1; i<=simul.get_num_cell();i++){
        B_p1[i] = dV[i]/(exp(dV[i])-1.0);
        B_p2[i] = B_p1[i]*exp(dV[i]);
    }
}
