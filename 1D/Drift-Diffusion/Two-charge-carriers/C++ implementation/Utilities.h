#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include "parameters.h"
#include <fstream>
#include <iomanip>
#include <iostream>

class Utilities
{
public:
    Utilities();

    std::vector<double> linear_mix(Parameters &params, std::vector<double> &new_values, std::vector<double> &old_values);
    void write_details(Parameters &params, double Va, std::vector<double> &V, std::vector<double> &p,  std::vector<double> &n, std::vector<double> &J_total, std::vector<double> &Un, std::vector<double> &PhotogenRate, std::vector<double> &R_Langevin);
    void write_JV(Parameters &params, std::ofstream &JV, double iter, double Va, std::vector<double> &J_total);

};

#endif // UTILITIES_H
