#include "optimization.h"
#include <random>


//constructor
Optim::Gradient_Descent::Gradient_Descent(Parameters &params): Params(params)  //NOTE: references must be initialized. //create a reference to the params structure so it can be used by all the member functions
{
    num_steps = 10;  //number of steps to take for each parameter--> determines the fine-ness of the optimization
    num_restarts = 3; //number of GD runs to do (each time starting with random numbers (within specified ranges) for the fitting variables. Is needed to avoid getting a local min.

    sign = 1;    //determines the sign for the parameter adjustment. By default start with +1.
    iter = 1;
    inner_iter=1;  //this is iter count for the steps at each order of magnitude vlue, i.e. when devide step size /10, this iter count refreshes

    overshoot = false;  //this becomes true when have overshot a local min--> needed to properly go back to the min

    best_lsqr_diff = 1e200;  //very large number to start

    best_vars.resize(Params.vars.size());

    //prepare for optimization:
    //read in experimental curve into vectors--> since it never changes, just do once

    exp_JV.open(params.exp_data_file_name);
    //check if file was opened
    if (!exp_JV) {
        std::cerr << "Unable to open file " << params.exp_data_file_name <<"\n";
        exit(1);   // call system to stop
    }

    double temp_V, temp_J;   //for input until end of file
    while (exp_JV >> temp_V >> temp_J) {  //there are 2 entries / line
        V_vector.push_back(temp_V);
        J_vector_exp.push_back(temp_J);
    }

}

void Optim::Gradient_Descent::run_GD()
{
    for (int GD_runs = 1; GD_runs < num_restarts; GD_runs++) {

        for (int var_index = 0; var_index < Params.vars.size(); var_index++) {   //for each variable that are optimizing

            //choose starting value randomly (need restarts to ensure that are not stuck in a local min)
            *Params.vars[var_index] = rand_num(Params.vars_min[var_index], Params.vars_max[var_index]);

            //reset iter counts
            iter = 1;
            inner_iter = 1;

            //bool min_found = false;  //used to decide whether to continue iterations for each variable
            bool fine_tuning = false; //to identify whether have already  did step/10 for fine tuning

            double step = (Params.vars_max[var_index]-Params.vars_min[var_index])/num_steps;  //step for the gradient descent
            //NOTE: there's risk of getting negative numbers when do this way, so need if's to protect from that

            do {

                if (overshoot) {
                    //don't over-write old_lsqr_diff--> let it be the prev. value b/c will return to that
                    old_lsqr_diff = old_lsqr_diff;
                    overshoot = false;  //turn it off, after corrected for overshoot
                } else {
                    old_lsqr_diff = lsqr_diff;  //save old value for figuring out if are moving towards or away from local minimum
                }

                //-----------------------------------------------------------------------------
                // Calculate the least squares difference between experimental and theory curve

                lsqr_diff = cost_function();  //this RUNS the DD simulation and computes the lsqr_diff

                //for test:
                std::cout << lsqr_diff << " at iteration " << iter << std::endl;
                std::cout << "variable # " << var_index << " value " << *Params.vars[var_index] << "sign " << sign << std::endl << std::endl;

                //std::cout << "inner iter " << inner_iter << std::endl;
                std::cout << "old diff " << old_lsqr_diff << std::endl;
                std::cout << "new lsqr diff " << lsqr_diff << std::endl << std::endl;

                //-----------------------------------------------------------------------------
                //Adjust the parameters

                //at 1st iter, always  start by going to in the direction determined by "sign". ASSUME THAT the selected initial guess value is roughly in the middle of the min, max rangee, since this is logical.
                if (inner_iter ==1) {
                    *Params.vars[var_index] += sign*step; //use dereference (*) to change the value that are pointing to
                }

                //at 2nd iter, need to decide which direction to go (i.e change parameter higher or lower)
                else if (inner_iter ==2) {
                    if (lsqr_diff < old_lsqr_diff) { // if less than, then means that moving to in the "sign" direction is reducing lsqr, so continue doing so
                        *Params.vars[var_index] += sign*step;
                    }
                    else{  //need to move to the left, past the initial guess, and 1 more step to left
                        *Params.vars[var_index] -= sign*2*step;  //2* b/c need to move back 1 and then 1 more to left
                        std::cout << "are movign to the left " <<  *Params.vars[var_index] << std::endl << std::endl;

                        //NOTE: IF I HAVE PARAMETERS WHICH ARE SUPPOSED TO BE NEGATIVE, THIS WILL NEED TO BE MODIFIED!!!
                        if (*Params.vars[var_index] < 0) {//took too large step, so reduce step size
                            std::cout << "value bacame negative, lower step size " << std::endl;
                            *Params.vars[var_index] += sign*step; //move back 1 step
                            step = step/10.;
                            *Params.vars[var_index] -= sign*step;
                        }
                        sign = -1*sign;  //NEED TO CHANGE THE SIGN, since now need to move the other way

                    }
                }

                //for all other iters
                else {
                    if (lsqr_diff < old_lsqr_diff) {
                        *Params.vars[var_index] += sign*step;

                        if (*Params.vars[var_index] < 0) {//took too large step, so reduce step size
                            std::cout << "value bacame negative, lower step size " << std::endl;
                            *Params.vars[var_index] -= sign*step; //other way 1 step
                            step = step/10.;
                            *Params.vars[var_index] += sign*step;
                        }
                        //CHECK for negative value again
                        if (*Params.vars[var_index] < 0) {
                            std::cout << "value AGAIN bacame negative. Will EXIT " << std::endl;
                            exit(1);
                        }
                    }

                    else if (lsqr_diff > old_lsqr_diff) {
                        *Params.vars[var_index] -= sign*step;
                        //if lsqr diff is greater in this iter than in previous iter, then means have already passed the minimum, so GO BACK TO previous value of
                        //the parameters, which is the local min value, and start optimizing other parameters

                        if(fine_tuning) {
                            lsqr_diff = old_lsqr_diff;  //b/c I step back 1 from where overshot...
                            std::cout << "Local min for variable # " << var_index << " = " << *Params.vars[var_index] << "has been found." << std::endl;
                            std::cout << "Cost at local min = " << lsqr_diff << std::endl; //=old_lsqr_diff b/c
                            break;           //come out of the while loop, if are fine tuning and lsqr_diff begins increasing, means that have found the local min
                        }
                        overshoot = true;  //to not rewrite the old_lsqr value,since taking 1 step back, using this overshoot condition

                        std::cout << "The approx. value corresponding to minimum least squares difference is " << *Params.vars[var_index] << std::endl;;
                        std::cout << "Will now fine tune the value " << std::endl;

                        //-----------------------------------------------------------------------------
                        //now take the step and divide by 10 again, for the future iterations..., will use the smaller step to fine tune the value

                        step = step/10.;
                        fine_tuning = true;
                        sign = 1;  //reset the sign to +1 for the new search, in the finer range
                        inner_iter = 1;  //need to reset inner_iter to 1 so can determine which direction to go.... again
                    }
                }

                iter++;
                inner_iter++;
                std::cout << "step size " << step << std::endl;

            } while (iter < Params.optim_max_iter);  //continue the iterations for a single parameter, while the lsqr_diff is decreasing, and are below max iter #

        }

        if (lsqr_diff < best_lsqr_diff) {
            best_lsqr_diff = lsqr_diff;
            std::cout << "A better parameter set has been found with cost = " << lsqr_diff << std::endl;
            std::cout << "Best parameter values are " << std::endl;

            for (int i = 0; i < Params.vars.size(); i++) {
                best_vars[i] = *Params.vars[i];
                std::cout << best_vars[i] << std::endl;
            }
        }
    }
}

//Particle struct constructor
Optim::Particle_swarm::Particle::Particle(int n_vars, Parameters &params)
{
    particle_params = params;   //make particle's parameters = to overall params. NOTE: this makes a copy,
    velocity.resize(n_vars);
    position.resize(n_vars);
    best_position.resize(n_vars);
    particle_vars.resize(n_vars);

    //initialize the POINTERS to particle_params for the vars that will need changing
    //FOR NOW THE  ONLY WAY I SEE TO DO IT IS HARDCODE THE PARAMETERS THAT NEED CHANGING...

    particle_vars[0] = &particle_params.Photogen_scaling;
    particle_vars[1] = &particle_params.n_mob_active;

}

//particle swarm constructor
Optim::Particle_swarm::Particle_swarm(Parameters &params) : Params(params)  //NOTE: references must be initialized. //create a reference to the params structure so it can be used by all the member functions
{

    //prepare for optimization:
    //read in experimental curve into vectors--> since it never changes, just do once

    exp_JV.open(params.exp_data_file_name);
    //check if file was opened
    if (!exp_JV) {
        std::cerr << "Unable to open file " << params.exp_data_file_name <<"\n";
        exit(1);   // call system to stop
    }

    double temp_V, temp_J;   //for input until end of file
    while (exp_JV >> temp_V >> temp_J) {  //there are 2 entries / line
        V_vector.push_back(temp_V);
        J_vector_exp.push_back(temp_J);
    }

    // Parameters of PSO (these might later need to be moved to the input file and params structure)
    //MOVE THE PARAMETERS TO AN INPUT FILE LATER.....
    PSO_max_iters = 100;        // Max # of iterations for PSO

    n_particles = 10;          // Population (Swarm) Size
    n_vars = Params.vars.size();   //number of variables that are adjusting

    //Clerc-Kennedy Constriction
    kappa = 1;
    phi1 = 2.05;
    phi2 = 2.05;
    phi = phi1 + phi2;
    chi = 2*kappa/abs(2-phi - sqrt(phi*phi - 4*phi));

    /*
    //PSO coefficients (WITHOUT Clerc-Kennedy constriction)
    w = 1;              // Inertia coefficient
    wdamp = 0.99;       // Damping for Inertia coefficient
    c1 = 2;             // Personal acceleration coefficient
    c2 = 2;             // Social (global) acceleration coefficient
    */

    //PSO coefficients WITH Clerc-Kennedy constriction
    w = chi;              // Inertia coefficient
    wdamp = 1;           // Damping for Inertia coefficient
    c1 = chi*phi1;       // Personal acceleration coefficient
    c2 = chi*phi2;             // Social (global) acceleration coefficient


    //Limit the velocity (for particle movements in parameter space)
    max_vel.resize(n_vars);
    min_vel.resize(n_vars);
    for (int i = 0; i < n_vars; i++) {
        max_vel[i] = 0.2*(Params.vars_max[i]-Params.vars_min[i]);  //the max and min velocities are dependent on the max and min values of variables
        min_vel[i] = -max_vel[i];
    }

    //Initialize global best
    global_best_cost = 1e200;    //since are minimizing, set global best to be very high

    //---------------------------------------------------------------------------------------
    //Initialization
    for (int i = 1; i <= n_particles; i++) {
       Particle *particle = new Particle(n_vars, Params);

       for (int dim = 0; dim < n_vars; dim++) {
            particle->position[dim] = rand_num(Params.vars_min[dim], Params.vars_max[dim]); //use uniformly distributed random number btw min and max value of the variable's range
            particle->velocity[dim] = 0;    


           //NOTE: for each particle, need to use a different parameter set which corresponds to that particle!
           //for now do like this--> but later need to make more elegant

            *particle->particle_vars[dim] = particle->position[dim];  //change the vars values (which are addresses of the particle_params which need changing)
            //std::cout << "particle position " << particle->position[dim] << std::endl;

       }

       particle->cost = cost_function(particle->particle_params); //run DD to find the cost function for the current particle (parameter set)

       //update the personal best
       particle->best_position = particle->position;
       particle->best_cost = particle->cost;

       //update global best
       if (particle->best_cost < global_best_cost) {
           global_best_cost = particle->best_cost;
           global_best_position = particle->best_position;
       }

       particles.push_back(particle);
    }
}



void Optim::Particle_swarm::run_PSO()
{
    //Main loop of PSO
    int iter = 0;

    do {
        iter++;

        for (int i = 0; i < n_particles; i++) {  //from 0 b/c of indexing
            Particle particle = *particles[i]; //need to use this to dereference the pointer and be able to access data memebers

            //Update Velocity
            for (int dim = 0; dim < n_vars; dim++) { //update needed for each dimension in cost space

                particle.velocity[dim] = w*particle.velocity[dim]
                        + c1*rand_num(0.0, 1.0)*(particle.best_position[dim] - particle.position[dim])
                        + c2*rand_num(0.0, 1.0)*(global_best_position[dim] - particle.position[dim]);

                //Update Position
                particle.position[dim] += particle.velocity[dim];

                //Apply Velocity Limits
                particle.velocity[dim] = std::max(particle.velocity[dim], min_vel[dim]);
                particle.velocity[dim] = std::min(particle.velocity[dim], max_vel[dim]);

                //Apply Lower and Upper Bound Limits (a clamping mechanism)
                particle.position[dim] = std::max(particle.position[dim], Params.vars_min[dim]); //if particle position is lower than VarMin, then the position becomes VarMin
                particle.position[dim] = std::min(particle.position[dim], Params.vars_max[dim]);


                //===================================================================================
                //Evaluation
                //NOTE: for each particle, need to use a different parameter set which corresponds to that particle!
                //for now do like this--> but later need to make more elegant


                *particle.particle_vars[dim] = particle.position[dim];  //change the vars values (which are addresses of the particle_params which need changing)
                particle.cost = cost_function(particle.particle_params);//call run_DD here with the current particle's parameters....;

                std::cout << "particle position " << dim << particle.position[dim] << std::endl;

            }

            //THIS NEEDS TO BE IMPROVED LATER TO MAKE IT MORE FLEXIBLE/GENERAL!

            //===================================================================================

            // Update Personal Best
            if (particle.cost < particle.best_cost) {
                 particle.best_position = particle.position;
                 particle.best_cost = particle.cost;

                 //Update Global Best
                 if (particle.best_cost < global_best_cost) {
                    global_best_cost = particle.best_cost;
                    global_best_position = particle.position;
                 }

            }
        }
        // Store the Best Cost Value at every iteration
        global_best_costs.push_back(global_best_cost);

        std::cout << "Iteration " << iter << ": Best Cost = " << global_best_costs[iter-1]<< std::endl;  //use iter -1 b/c I'm counting iters from 1, but vectors index from 0
        std::cout << "Photogenrate value " << global_best_position[0] << std::endl;
        std::cout << "n_mob value " << global_best_position[1] << std::endl;

        //Damp Inertial  Coefficient in each iteration (note: only used when not using Clerc-Kennedy restriction)
        w = w * wdamp;

    } while (global_best_cost > Params.fit_tolerance && iter < PSO_max_iters);
}

//--------------------------------------------------------------------------------------------------------------------------
//generates a uniformly  distributed random number in the range [a,b]
double Optim::Particle_swarm::rand_num(double a, double b)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::default_random_engine generator(seed);   //NOTE: need static, so is  called only on 1st call (otherwise, since random #'s are pseudo random, will get same # on each call
    std::uniform_real_distribution<double> distribution(a, b);

    return distribution(generator);
}

double Optim::Gradient_Descent::rand_num(double a, double b)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::default_random_engine generator(seed);   //NOTE: need static, so is  called only on 1st call (otherwise, since random #'s are pseudo random, will get same # on each call
    std::uniform_real_distribution<double> distribution(a, b);

    return distribution(generator);
}

//-----------------------------------------------------------------------------------------------------------------------------------

//Note: this is not the most elegant way, but for now we use a get_exp_data function definition within each optimization method's class

void Optim::Gradient_Descent::get_exp_data()
{
    //prepare for optimization:
    //read in experimental curve into vectors--> since it never changes, just do once

    exp_JV.open(Params.exp_data_file_name);
    //check if file was opened
    if (!exp_JV) {
        std::cerr << "Unable to open file " << Params.exp_data_file_name <<"\n";
        exit(1);   // call system to stop
    }

    double temp_V, temp_J;   //for input until end of file
    while (exp_JV >> temp_V >> temp_J) {  //there are 2 entries / line
        V_vector.push_back(temp_V);
        J_vector_exp.push_back(temp_J);
    }
}

void Optim::Particle_swarm::get_exp_data()
{
    //prepare for optimization:
    //read in experimental curve into vectors--> since it never changes, just do once

    exp_JV.open(Params.exp_data_file_name);
    //check if file was opened
    if (!exp_JV) {
        std::cerr << "Unable to open file " << Params.exp_data_file_name <<"\n";
        exit(1);   // call system to stop
    }

    double temp_V, temp_J;   //for input until end of file
    while (exp_JV >> temp_V >> temp_J) {  //there are 2 entries / line
        V_vector.push_back(temp_V);
        J_vector_exp.push_back(temp_J);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------

double Optim::Gradient_Descent::cost_function()
{
   J_vector_model = run_DD(Params);
   lsqr_diff = 0;  //rezero
   //calculate least squares
   //just compare the J values (V values for experiment and theory should be the same)
   for (int i = 0; i < J_vector_model.size(); i++) {
       lsqr_diff += (J_vector_model[i] - J_vector_exp[i])*(J_vector_model[i] - J_vector_exp[i]);
   }

   return lsqr_diff;
}


double Optim::Particle_swarm::cost_function(Parameters &particle_params)
{
    //NOTE: need to run DD based on the particle positions..., so Params should be different for each particle....
   J_vector_model = run_DD(particle_params);
   lsqr_diff = 0;  //rezero
   //calculate least squares
   //just compare the J values (V values for experiment and theory should be the same)
   for (int i =0; i < J_vector_model.size(); i++) {
       lsqr_diff += (J_vector_model[i] - J_vector_exp[i])*(J_vector_model[i] - J_vector_exp[i]);
   }

   return lsqr_diff;
}
