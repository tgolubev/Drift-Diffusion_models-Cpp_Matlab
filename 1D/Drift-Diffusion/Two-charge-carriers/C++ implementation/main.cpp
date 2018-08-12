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
#include "optimization.h"

int main()
{

    Parameters params;    //params is struct storing all parameters
    params.Initialize();  //reads parameters from file (will set the starting parameter set for PSO)


    if (params.auto_fit == true) {

        if (params.optim_method == 1) {
            Optim::Gradient_Descent GD(params);
            GD.run_GD();
        } else if (params.optim_method ==2) {
            Optim::Particle_swarm PSO(params);  //creating an object of Particle_swarm struct which is member of Optimization class
            PSO.run_PSO();
        } else {
            std::cout << "Invalid optimization method" << std::endl;
            exit(1);
        }

    } else {
        run_DD(params);  //run once
    }


    return 0;
}

