#include "poisson.h"
#include <iostream>

//---------------Setup AV diagonals (Poisson solve)---------------------------------------------------------------
//this is a in tridiag_solver
void Poisson::set_main_diag(const std::vector<double> &epsilon){
    int num_elements = main_diag.size() - 1; //num_elements corresponds to the particular diagonal. -1 b/c we don't use the 0th element.

    for (int i=1; i<= num_elements;i++){
        main_diag[i] = -2.*epsilon[i];
    }
}

//this is b in tridiag_solver
void Poisson::set_upper_diag(const std::vector<double> &epsilon){
    int num_elements = upper_diag.size() -1;  //-1 since vectors fill from 0, but I'm indexing from 1
    //NOTE: num_elements here, is actually # of elements in the off-diagonal

    for (int i = 1; i<=num_elements; i++){
        upper_diag[i] = epsilon[i];
    }
}

//this is c in tridiag_solver
void Poisson::set_lower_diag(const std::vector<double> &epsilon){
    int num_elements = lower_diag.size() -1;

    for (int i = 1; i<=num_elements; i++){
        lower_diag[i] = epsilon[i];
    }
}

//------------------------------------------------------------------------------------
void Poisson::set_rhs(const std::vector<double> &epsilon, const std::vector<double> &n, const std::vector<double> &p, double V_leftBC, double V_rightBC){
    for(int i = 1;i<=rhs.size()-1; i++){
        rhs[i] = CV*(n[i] - p[i]);
        //std::cout << "bV " << bV[i] << std::endl;
    }
    rhs[1] -= epsilon[0]*V_leftBC;
    rhs[rhs.size()-1] -= epsilon[rhs.size()]*V_rightBC;
}
