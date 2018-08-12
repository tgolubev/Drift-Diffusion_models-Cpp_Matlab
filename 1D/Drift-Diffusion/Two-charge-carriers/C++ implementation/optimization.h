#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>   //allows to use fill and min
#include <fstream>
#include <chrono>
#include <string>
#include <time.h>
#include <fstream>
#include <string>

#include "run_DD.h"
#include "parameters.h"

class Optimization
{
public:
    Optimization::Optimization(Parameters &params);
    void gradient_descent(Parameters &params);
    void particle_swarm(Parameters &params);
private:


    std::vector<double> V_vector;
    std::vector<double> J_vector_exp;
    double lsqr_diff, old_lsqr_diff;
    std::vector<double> J_vector_model;

    std::ifstream exp_JV;


};



#endif // OPTIMIZATION_H
