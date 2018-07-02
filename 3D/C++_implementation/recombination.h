#ifndef RECOMBINATION_H
#define RECOMBINATION_H

#include<vector>
#include "constants.h"
#include "parameters.h"

class Recombo{
public:
    //!Constructor sets up the recombination member parameters using the input from \param params
    Recombo(const Parameters &params);

    //!Computes bimolecular Langevin recombination rate.
    //! This depends on the electron \param n and hole \param p density.
    std::vector<double> ComputeR_Langevin(const Parameters &params, const std::vector<double> &n, const std::vector<double> &p);

private:
    double k_rec; //!bimolecular recombination coefficient
    double E_trap; //!Trap level energy
    double n1;
    double p1;
    std::vector<double> R_Langevin;  //!store the Langevin recombination rate
};


#endif // RECOMBINATION_H
