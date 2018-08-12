#include "optimization.h"


//constructor
Optim::Gradient_Descent::Gradient_Descent(Parameters &params): Params(params)  //NOTE: references must be initialized. //create a reference to the params structure so it can be used by all the member functions
{
    num_steps = 10;  //number of steps to take for each parameter--> determines the fine-ness of the optimization
    sign = 1;    //determines the sign for the parameter adjustment. By default start with +1.
    iter = 1;
    inner_iter=1;  //this is iter count for the steps at each order of magnitude vlue, i.e. when devide step size /10, this iter count refreshes
    Photogen_scaling_step = (params.Photogen_scaling_max-params.Photogen_scaling_min)/num_steps;
    overshoot = false;  //this becomes true when have overshot a local min--> needed to properly go back to the min

}

void Optim::Gradient_Descent::run_GD()
{
    do {

        if (overshoot)
            //don't over-write old_lsqr_diff--> let it be the prev. value b/c will return to that
            old_lsqr_diff = old_lsqr_diff;
        else
            old_lsqr_diff = lsqr_diff;  //save old value for figuring out if are moving towards or away from local minimum

        lsqr_diff = 0;  //rezero

        //-----------------------------------------------------------------------------
        // Calculate the least squares difference between experimental and theory curve

        lsqr_diff = cost_function();  //this RUNS the DD simulation and computes the lsqr_diff

        //for test:
        std::cout << lsqr_diff << " at iteration " << iter << std::endl;
        std::cout << "photogenscalling" << Params.Photogen_scaling << "sign " << sign << std::endl;

        //-----------------------------------------------------------------------------
        //Adjust the parameters

        //at 1st iter, always  start by going to in the direction determined by "sign". ASSUME THAT the selected initial guess value is roughly in the middle of the min, max rangee, since this is logical.
        if (inner_iter ==1) {
            Params.Photogen_scaling += sign*Photogen_scaling_step;
        }

        //at 2nd iter, need to decide which direction to go (i.e change parameter higher or lower)
        else if (inner_iter ==2) {
            if (lsqr_diff < old_lsqr_diff) { // if less than, then means that moving to in the "sign" direction is reducing lsqr, so continue doing so
                Params.Photogen_scaling += sign*Photogen_scaling_step;
            }
            else{  //need to move to the left, past the initial guess, and 1 more step to left
                Params.Photogen_scaling -= sign*2*Photogen_scaling_step;  //2* b/c need to move back 1 and then 1 more to left
                sign = -1*sign;  //NEED TO CHANGE THE SIGN, since now need to move the other way
            }
        }

        //for all other iters
        else {
            if (lsqr_diff < old_lsqr_diff)
                Params.Photogen_scaling += sign*Photogen_scaling_step;     //we will optimize Photogeneration-scaling parameter for the test
            else if (lsqr_diff > old_lsqr_diff) {
                Params.Photogen_scaling -= sign*Photogen_scaling_step;
                //if lsqr diff is greater in this iter than in previous iter, then means have already passed the minimum, so GO BACK TO previous value of
                //the parameters, which is the local min value, and start optimizing other parameters

                overshoot = true;  //to not rewrite the old_lsqr value,since taking 1 step back, using this overshoot condition

                std::cout << "The value corresponding to minimum least squares difference is " << Params.Photogen_scaling << std::endl;;
                std::cout << "Will now fine tune the value " << std::endl;

                //-----------------------------------------------------------------------------
                //now take the step and divide by 10 again, for the future iterations..., will use the smaller step to fine tune the value

                Photogen_scaling_step = Photogen_scaling_step/10.;
                std::cout << Photogen_scaling_step << std::endl;
                sign = 1;  //reset the sign to +1 for the new search, in the finer range
                inner_iter = 1;  //need to reset inner_iter to 1 so can determine which direction to go.... again
            }
        }

        iter++;
        inner_iter++;

    } while (lsqr_diff > Params.fit_tolerance && iter < Params.optim_max_iter);
}

//particle swarm constructor
Optim::Particle_swarm::Particle_swarm(Parameters &params) : Params(params)  //NOTE: references must be initialized. //create a reference to the params structure so it can be used by all the member functions
{

    // Parameters of PSO (these might later need to be moved to the input file and params structure)
    max_iters = 1000;        // Max # of iterations for PSO

    n_particles = 50;          // Population (Swarm) Size
    n_vars = 1;   //number of variables that are adjusting

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

    //setup the vector with min and max values  for positions (the variables)
    var_min.resize(n_vars);
    var_max.resize(n_vars);
    var_min[0] = params.Photogen_scaling_min;
    var_max[0] = params.Photogen_scaling_max;

    //WILL NEED TO ADD OTHER LIMITS TO PARAMETERS WHEN OPTIMIZE FOR MORE PARAMETERS

    //Limit the velocity (for particle movements in parameter space)
    max_vel.resize(n_vars);
    min_vel.resize(n_vars);
    for (int i = 0; i < n_vars; i++) {
        max_vel[i] = 0.2*(var_max[i]-var_min[i]);  //the max and min velocities are dependent on the max and min values of variables
        min_vel[i] = -max_vel[i];
    }


    //Initialize global best
    global_best_cost = 1e200;    //since are minimizing, set global best to be very high


    //---------------------------------------------------------------------------------------
    //Initialization
    for (int i = 1; i <= n_particles; i++) {
       Particle *particle = new Particle();
       particle->velocity.resize(n_vars);

       for (int dim = 0; dim <= n_vars; dim++) {
            particle->position[dim] = rand();//CALL RANDOM NUMBER GENERATOR
            particle->velocity[dim] = 0;
       }

       particle->cost = cost_function(); //run DD to find the cost function for the current particle (parameter set)

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
};



void Optim::Particle_swarm::run_PSO()
{
    //Main loop of PSO

    for (int it = 1; it <= max_iters; it++) {
        for (int i = 0; i < n_particles; i++) {  //from 0 b/c of indexing
            Particle particle = *particles[i]; //need to use this to dereference the pointer and be able to access data memebers

            //Update Velocity
            for (int dim = 0; dim < n_vars; dim++) { //update needed for each dimension in cost space

                particle.velocity[dim] = w*particle.velocity[dim]
                        + c1*rand()*(particle.best_position[dim] - particle.position[dim])
                        + c2*rand()*(global_best_position[dim] - particle.position[dim]);

                //Update Position
                particle.position[dim] += particle.velocity[dim];

                //Apply Velocity Limits
                particle.velocity[dim] = std::max(particle.velocity[dim], min_vel[dim]);
                particle.velocity[dim] = std::min(particle.velocity[dim], max_vel[dim]);

                //Apply Lower and Upper Bound Limits (a clamping mechanism)
                particle.position[dim] = std::max(particle.position[dim], var_min[dim]); //if particle position is lower than VarMin, then the position becomes VarMin
                particle.position[dim] = std::min(particle.position[dim], var_max[dim]);
            }

            //Evaluation
            particle.cost = cost_function();//call run_DD here with the current particle's parameters....;

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

        std::cout << "Iteration " << ": Best Cost = " << global_best_costs[it];

        //Damp Inertial  Coefficient in each iteration (note: only used when not using Clerc-Kennedy restriction)
        w = w * wdamp;
    }
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

   //calculate least squares
   //just compare the J values (V values for experiment and theory should be the same)
   for (int i =0; i < J_vector_model.size(); i++) {
       lsqr_diff += (J_vector_model[i] - J_vector_exp[i])*(J_vector_model[i] - J_vector_exp[i]);
   }

   return lsqr_diff;
}


double Optim::Particle_swarm::cost_function()
{
   J_vector_model = run_DD(Params);

   //calculate least squares
   //just compare the J values (V values for experiment and theory should be the same)
   for (int i =0; i < J_vector_model.size(); i++) {
       lsqr_diff += (J_vector_model[i] - J_vector_exp[i])*(J_vector_model[i] - J_vector_exp[i]);
   }

   return lsqr_diff;
}
