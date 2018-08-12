#include "optimization.h"

//constructor
Optimization::Optimization(Parameters &params)
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



}

void Optimization::gradient_descent(Parameters &params)
{
    int num_steps = 10;  //number of steps to take for each parameter--> determines the fine-ness of the optimization
    int sign = 1;    //determines the sign for the parameter adjustment. By default start with +1.


    int iter = 1;
    int inner_iter=1;  //this is iter count for the steps at each order of magnitude vlue, i.e. when devide step size /10, this iter count refreshes
    double Photogen_scaling_step = (params.Photogen_scaling_max-params.Photogen_scaling_min)/num_steps;
    bool overshoot = false;  //this becomes true when have overshot a local min--> needed to properly go back to the min

    do {

        J_vector_model = run_DD(params);  //in future this here can be distributed among many processors....

        //-----------------------------------------------------------------------------
        // Calculate the least squares difference between experimental and theory curve

        if (overshoot)
            //don't over-write old_lsqr_diff--> let it be the prev. value b/c will return to that
            old_lsqr_diff = old_lsqr_diff;
        else
            old_lsqr_diff = lsqr_diff;  //save old value for figuring out if are moving towards or away from local minimum

        lsqr_diff = 0;  //rezero
        //just compare the J values (V values for experiment and theory should be the same)
        for (int i =0; i < J_vector_model.size(); i++) {
            lsqr_diff += (J_vector_model[i] - J_vector_exp[i])*(J_vector_model[i] - J_vector_exp[i]);
        }

        //for test:
        std::cout << lsqr_diff << " at iteration " << iter << std::endl;
        std::cout << "photogenscalling" << params.Photogen_scaling << "sign " << sign << std::endl;

        //-----------------------------------------------------------------------------
        //Adjust the parameters

        //at 1st iter, always  start by going to in the direction determined by "sign". ASSUME THAT the selected initial guess value is roughly in the middle of the min, max rangee, since this is logical.
        if (inner_iter ==1) {
            params.Photogen_scaling += sign*Photogen_scaling_step;
        }

        //at 2nd iter, need to decide which direction to go (i.e change parameter higher or lower)
        else if (inner_iter ==2) {
            if (lsqr_diff < old_lsqr_diff) { // if less than, then means that moving to in the "sign" direction is reducing lsqr, so continue doing so
                params.Photogen_scaling += sign*Photogen_scaling_step;
            }
            else{  //need to move to the left, past the initial guess, and 1 more step to left
                params.Photogen_scaling -= sign*2*Photogen_scaling_step;  //2* b/c need to move back 1 and then 1 more to left
                sign = -1*sign;  //NEED TO CHANGE THE SIGN, since now need to move the other way
            }
        }

        //for all other iters
        else {
            if (lsqr_diff < old_lsqr_diff)
                params.Photogen_scaling += sign*Photogen_scaling_step;     //we will optimize Photogeneration-scaling parameter for the test
            else if (lsqr_diff > old_lsqr_diff) {
                params.Photogen_scaling -= sign*Photogen_scaling_step;
                //if lsqr diff is greater in this iter than in previous iter, then means have already passed the minimum, so GO BACK TO previous value of
                //the parameters, which is the local min value, and start optimizing other parameters

                overshoot = true;  //to not rewrite the old_lsqr value,since taking 1 step back, using this overshoot condition

                std::cout << "The value corresponding to minimum least squares difference is " << params.Photogen_scaling << std::endl;;
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

    } while (lsqr_diff > params.fit_tolerance && iter < params.optim_max_iter);

}

void Optimization::particle_swarm(Parameters &params)
{




}
