#include "parameters.h"
#include <string>


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

    try{
        std::string comment;  //to "eat" the comments. will only work is comment has no spaces btw words
        parameters >> comment;  //header line
        parameters >> L >> comment;
        isPositive(L, comment);
        parameters >> num_cell >> comment;
        isPositive(num_cell, comment);
        parameters >> N_LUMO >> comment;  //we will just ignore the comments
        isPositive(N_LUMO,comment);
        parameters >> N_HOMO >> comment;
        isPositive(N_HOMO,comment);
        parameters >> Photogen_scaling >> comment;
        isPositive(Photogen_scaling,comment);
        parameters >> phi_a >> comment;
        parameters >> phi_c >> comment;
        parameters >> eps_active >> comment;
        isPositive(eps_active,comment);
        parameters >> p_mob_active >> comment;
        isPositive(p_mob_active,comment);
        parameters >> n_mob_active >> comment;
        isPositive(n_mob_active,comment);
        parameters >> mobil >> comment;
        isPositive(mobil,comment);
        parameters >> E_gap >> comment;
        isPositive(E_gap,comment);
        parameters >> active_CB >> comment;  //this should be negative
        isNegative(active_CB, comment);
        parameters >> active_VB >> comment;  //should be negative
        isNegative(active_VB, comment);
        parameters >> WF_anode >> comment;
        isPositive(WF_anode ,comment);
        parameters >> WF_cathode >> comment;
        isPositive(WF_cathode,comment);
        parameters >> k_rec_input >> comment;
        isPositive(k_rec_input,comment);
        parameters >> dx >> comment;
        isPositive(dx ,comment);
        parameters >> Va_min >> comment;
        parameters >> Va_max >> comment;
        parameters >> increment >> comment;
        isPositive(increment,comment);
        parameters >> w_eq >> comment;
        isPositive(w_eq,comment);
        parameters >> w_i >> comment;
        isPositive(w_i,comment);
        parameters >> tolerance_i >> comment;
        isPositive(tolerance_i,comment);
        parameters >> w_reduce_factor >> comment;
        isPositive(w_reduce_factor,comment);
        parameters >> tol_relax_factor >> comment;
        isPositive(tol_relax_factor,comment);
        parameters.close();
        N = N_HOMO;     //scaling factor helps CV be on order of 1

    }
    catch(std::exception &e){
        std::cerr << e.what() << std::endl;
        exit(1);
    }


}

void Parameters::isPositive(double input, std::string comment)
{
    if(input <=0){
        std::cerr << "error: Non-positive input for " << comment << std::endl;
        std::cerr << "Input was read as " << input << std::endl;
        throw std::runtime_error("Invalid input. This input must be positive.");
    }
}

void Parameters::isPositive(int input, std::string comment)
{
    if(input <=0){
        std::cerr << "error: Non-positive input for " << comment << std::endl;
        std::cerr << "Input was read as " << input << std::endl;
        throw std::runtime_error("Invalid input. This input must be positive.");
    }
}

void Parameters::isNegative(double input, std::string comment)
{
    if(input >=0){
        std::cerr << "error: Non-negative input for " << comment << std::endl;
        std::cerr << "Input was read as " << input << std::endl;
        throw std::runtime_error("Invalid input. This input must be negative.");
    }
}

void Parameters::isNegative(int input, std::string comment)
{
    if(input >=0){
        std::cerr << "error: Non-negative input for " << comment << std::endl;
        std::cerr << "Input was read as " << input << std::endl;
        throw std::runtime_error("Invalid input. This input must be negative.");
    }

}
