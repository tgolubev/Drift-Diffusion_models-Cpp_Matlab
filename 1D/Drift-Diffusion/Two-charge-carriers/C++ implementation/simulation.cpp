#include "simulation.h"
#include "parameters.h"

#include<cmath>

Simulation::Simulation()      //this will read from file and setup all of the simulation parameters
{
  //for now just specify all parameters here...

    //Voltage sweep setup
    Va_min = -0.5;
    Va_max = 1.2;
    increment = 0.01;
    num_V = floor((Va_max-Va_min)/increment)+1;

    //Simulation parameters
    w_eq = 0.01;                     //for 1st equil convergence weighting factor
    w_i = 0.2;                       // %set up of weighting factor
    tolerance_eq = 100*5e-12;        //error tolerance for equil run
    tolerance_i = 5e-12;             //initial error tolerance for 1st non-equil Va (1s tVa with generation)d\

    w_reduce_factor = 2.;
    tol_relax_factor = 10.;

    //geometry parameters
    //Thicknesses (in m)
    L = 300.0e-9;             //device length in meters. Note: if change this and using a gen_rate file, need to update the number of elements in file to = num_cell-2
    dx = 1.e-9;
    num_cell = L/dx;            //number of cells

}
