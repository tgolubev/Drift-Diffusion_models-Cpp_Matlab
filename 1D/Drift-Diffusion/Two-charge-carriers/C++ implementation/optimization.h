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


    class Particle_swarm
    {
     public:
        Particle_swarm::Particle_swarm();
        void run_PSO(Parameters &params);

     private:
        std::vector<Particle*> particles;  //vector of pointers to particles

        double cost_function(Parameters &params);


        double global_best_cost;
        std::vector<double> global_best_position;

        //max and min velocity values (for PSO particle moves)
        std::vector<double> min_vel;
        std::vector<double> max_vel;

        //max and min  values for each of variables that are optimizinng
        std::vector<double> var_min;
        std::vector<double> var_max;

        int max_iters;        // Max # of iterations for PSO

        int n_particles;          // Population (Swarm) Size
        int n_vars;   //number of variables that are adjusting

        //Clerc-Kennedy Constriction
        double kappa;
        double phi1;
        double phi2;
        double phi;
        double chi;

        std::vector<double> global_best_costs;  //vector to store the best costs at each PSO iteration

        struct Particle
        {
            Particle::Particle();
            std::vector<double> position;
            std::vector<double> velocity;
            std::vector<double> cost;   //cost (aka objective) function which are minimizing, in our case will be least squares difference
            std::vector<double> best_position;
            double best_cost;

        };


    };



private:


    std::vector<double> V_vector;
    std::vector<double> J_vector_exp;
    double lsqr_diff, old_lsqr_diff;
    std::vector<double> J_vector_model;

    std::ifstream exp_JV;


};



#endif // OPTIMIZATION_H
