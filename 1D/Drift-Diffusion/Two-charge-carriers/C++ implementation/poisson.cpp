#include "poisson.h"
#include <iostream>

Poisson::Poisson(const Parameters &params)
{
    main_diag.resize(params.num_cell);
    upper_diag.resize(params.num_cell-1);
    lower_diag.resize(params.num_cell-1);
    rhs.resize(params.num_cell);
    epsilon.resize(params.num_cell +1);

    CV = params.N*params.dx*params.dx*q/(epsilon_0*Vt);
    std::fill(epsilon.begin(), epsilon.end(), params.eps_active);
} //constructor)


void Poisson::setup_matrix()  //Note: this is on purpose different than the setup_eqn used for Continuity eqn's, b/c I need to setup matrix only once
{
    set_main_diag();
    set_lower_diag();
    set_upper_diag();
}


//---------------Setup AV diagonals (Poisson solve)---------------------------------------------------------------
//this is a in tridiag_solver
void Poisson::set_main_diag()
{
    for (int i = 1; i < main_diag.size(); i++) {
        main_diag[i] = -2.*epsilon[i];
    }
}

//this is b in tridiag_solver
void Poisson::set_upper_diag()
{
    for (int i = 1; i < upper_diag.size(); i++) {
        upper_diag[i] = epsilon[i];
    }
}

//this is c in tridiag_solver
void Poisson::set_lower_diag()
{
    for (int i = 1; i < lower_diag.size(); i++) {
        lower_diag[i] = epsilon[i];
    }
}

//------------------------------------------------------------------------------------
void Poisson::set_rhs(const std::vector<double> &n, const std::vector<double> &p, double V_leftBC, double V_rightBC)
{
    for (int i = 1; i <= rhs.size()-1; i++) {
        rhs[i] = CV*(n[i] - p[i]);
        //std::cout << "bV " << bV[i] << std::endl;
    }
    rhs[1] -= epsilon[0]*V_leftBC;
    rhs[rhs.size()-1] -= epsilon[rhs.size()]*V_rightBC;
}
