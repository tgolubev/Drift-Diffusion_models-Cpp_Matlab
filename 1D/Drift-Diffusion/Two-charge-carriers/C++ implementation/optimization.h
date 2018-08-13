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

//This provides the namespace and class definitions

namespace Optim  //group all the optimization classes
{


    class Gradient_Descent
    {
    public:
        Gradient_Descent(Parameters &params);
        void run_GD();
        double cost_function();
        void get_exp_data();
        double rand_num(double a, double b);   //SHOULD later make this function defined somewhere like a new Utilities class for the code

    private:
        Parameters &Params; //reference to params so all member functions can use..

        int n_vars;   //number of variables that are adjusting

        int num_steps;  //number of steps to take for each parameter--> determines the fine-ness of the optimization
        int sign;    //determines the sign for the parameter adjustment. By default start with +1.


        int iter;
        int inner_iter;  //this is iter count for the steps at each order of magnitude vlue, i.e. when devide step size /10, this iter count refreshes
        double Photogen_scaling_step;
        bool overshoot;  //this becomes true when have overshot a local min--> needed to properly go back to the min
        int num_restarts;

        std::vector<double> best_vars;

        std::vector<double> V_vector;
        std::vector<double> J_vector_exp;
        double lsqr_diff, old_lsqr_diff, best_lsqr_diff;
        std::vector<double> J_vector_model;
        std::ifstream exp_JV;

    };


    //----------------------------------------------------------------------------------------------------
    class Particle_swarm
    {
    public:
        Particle_swarm(Parameters &params);
        void run_PSO();
        double cost_function(Parameters &particle_params);
        void get_exp_data();
        double rand_num(double a, double b);

    private:

        Parameters &Params; //reference to params so all member functions can use..

        double global_best_cost;
        std::vector<double> global_best_position;

        //max and min velocity values (for PSO particle moves)
        std::vector<double> min_vel;
        std::vector<double> max_vel;

        //max and min  values for each of variables that are optimizinng
        std::vector<double> var_min;
        std::vector<double> var_max;

        int PSO_max_iters;        // Max # of iterations for PSO

        int n_particles;          // Population (Swarm) Size
        int n_vars;   //number of variables that are adjusting

        //Clerc-Kennedy Constriction
        double kappa;
        double phi1;
        double phi2;
        double phi;
        double chi;

        double w, wdamp, c1, c2;

        std::vector<double> global_best_costs;  //vector to store the best costs at each PSO iteration

        struct Particle
        {
            Particle::Particle(int n_vars, Parameters &params);
            Parameters particle_params;   //NOTE: this is a COPY, not a reference, since want to change each particle's params individually.
            std::vector<double> position;
            std::vector<double> velocity;
            double cost;   //cost (aka objective) function which are minimizing, in our case will be least squares difference
            std::vector<double> best_position;
            double best_cost;

        };

        std::vector<Particle*> particles;  //NOTE: this declaration must follow after the struct Particle definition, otherwise, it doesn't know what a Paricle is!//vector of pointers to particles

        std::vector<double> V_vector;
        std::vector<double> J_vector_exp;
        double lsqr_diff, old_lsqr_diff;
        std::vector<double> J_vector_model;


        std::ifstream exp_JV;

    };


}

#endif // OPTIMIZATION_H
