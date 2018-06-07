#ifndef RECOMBINATION_H
#define RECOMBINATION_H

#include<vector>
#include "constants.h"
#include "parameters.h"


class Recombo{
public:
    Recombo(Parameters &params, double k_rec_input);  //constructor

    std::vector<double> ComputeR_Langevin(Parameters &params, const std::vector<double> &n, const std::vector<double> &p);
    //void set_k_rec(double k_rec_input){ k_rec = k_rec_input;}
private:
    double k_rec;
    double E_trap;
    double n1;
    double p1;
    std::vector<double> R_Langevin;
};


#endif // RECOMBINATION_H
