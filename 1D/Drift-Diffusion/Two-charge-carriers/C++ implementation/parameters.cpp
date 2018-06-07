#include "parameters.h"

#include<cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

Parameters::Parameters ()      //this will read from file and setup all of the Parameters parameters
{


}

void Parameters::Initialize()
{
    //input from file the parameters
    std::ifstream parameters;
    parameters.open("parameters.txt");
    //check if file was opened
    if (!parameters) {
        std::cerr << "Unable to open file parameters.txt";
        exit(1);   // call system to stop
    }

    std::string comment;  //to "eat" the comments. will only work is comment has no spaces btw words
    parameters >> comment;  //header line
    parameters >> L >> comment;
    parameters >> num_cell >> comment;
    parameters >> N_LUMO >> comment;  //we will just ignore the comments
    parameters >> N_HOMO >> comment;
    parameters >> Photogen_scaling >> comment;
    parameters >> phi_a >> comment;
    parameters >> phi_c >> comment;
    parameters >> eps_active >> comment;
    parameters >> p_mob_active >> comment;
    parameters >> n_mob_active >> comment;
    parameters >> mobil >> comment;
    parameters >> E_gap >> comment;
    parameters >> active_CB >> comment;
    parameters >> active_VB >> comment;
    parameters >> WF_anode >> comment;
    parameters >> WF_cathode >> comment;
    parameters >> k_rec_input >> comment;
    parameters >> dx >> comment;
    parameters >> Va_min >> comment;
    parameters >> Va_max >> comment;
    parameters >> increment >> comment;
    parameters >> w_eq >> comment;
    parameters >> w_i >> comment;
    parameters >> tolerance_i >> comment;
    parameters >> w_reduce_factor >> comment;
    parameters >> tol_relax_factor >> comment;
    parameters.close();
    N = N_HOMO;     //scaling factor helps CV be on order of 1


}
