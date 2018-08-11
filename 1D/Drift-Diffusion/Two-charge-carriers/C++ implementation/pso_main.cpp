/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Solving 1D Poisson + Drift Diffusion semiconductor eqns for a solar cell using
%                      Scharfetter-Gummel discretization

%  Automatic fitting by  PARTICLE SWARM OPTIMIZATION METHOD
%
%                         Written by Timofey Golubev
%
%     The code as is will calculate data for a JV curve
%     as well as carrier densities, current densities, and electric field
%     distributions of a generic solar cell made of an active layer and electrodes.
%     More equations for carrier recombination can be easily added.f
%
%     Photogeneration rate will be inputed from gen_rate.inp file
%     (i.e. the output of an optical model can be used) or an analytic expression
%     for photogeneration rate can be added to photogeneration.cpp. Generation rate file
%     should contain num_cell-2 number of entries in a single column, corresponding to
%     the the generation rate at each mesh point (except the endpoints).
%
%     The code can also be applied to non-illuminated devices by
%     setting photogeneration rate to 0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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

int main()
{

    Parameters params;    //params is struct storing all parameters
    params.Initialize();  //reads parameters from file (will set the starting parameter set for PSO)
    double lsqr_diff;
    std::vector<double> J_vector_model;

    //make the auto fitting an option in params--> if yes, do a loop, if no, then make no loop


    if (params.auto_fit == true) {

        //read in experimental curve into vectors--> since it never changes, just do once
        std::vector<double> V_vector;
        std::vector<double> J_vector_exp;
        std::ifstream exp_JV;
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

        //-----------------------------------------------------------------------------
        // do optimization
        int iter = 1;
        do {

            J_vector_model = run_DD(params);  //in future this here can be distributed among many processors....

            //-----------------------------------------------------------------------------
            // Calculate the least squares difference between experimental and theory curve

            //just compare the J values (V values for experiment and theory should be the same)
            for (int i =0; i < J_vector_model.size(); i++) {
                lsqr_diff += (J_vector_model[i] - J_vector_exp[i])*(J_vector_model[i] - J_vector_exp[i]);
            }

            //for test:
            std::cout << lsqr_diff << std::endl;

            //-----------------------------------------------------------------------------
            //Adjust the parameters


            iter++;

        } while (lsqr_diff > params.fit_tolerance && iter < params.optim_max_iter);





    } else {
        run_DD(params);  //run once
    }





    return 0;
}

