#include "continuity_n.h"

//-------------------------------Setup An diagonals (Continuity/drift-diffusion solve)-----------------------------
void Continuity_n::set_main_diag(const std::vector<double> &n_mob,const  std::vector<double> &B_n1,const  std::vector<double> &B_n2){
    int num_elements = main_diag.size() - 1;

    for (int i=1; i<= num_elements;i++){
        main_diag[i] = -(n_mob[i]*B_n1[i] + n_mob[i+1]*B_n2[i+1]);
    }
}

//this is b in tridiag_solver
void Continuity_n::set_upper_diag(const std::vector<double> &n_mob,const  std::vector<double> &B_n1){
    int num_elements = upper_diag.size() -1;  //-1 since vectors fill from 0, but I'm indexing from 1
    //NOTE: num_elements here, is actually # of elements in the off-diagonal

    for (int i = 1; i<=num_elements; i++){
        upper_diag[i] = n_mob[i+1]*B_n1[i+1];
    }
}

//this is c in tridiag_solver
void Continuity_n::set_lower_diag(const std::vector<double> &n_mob,const  std::vector<double> &B_n2){
    int num_elements = lower_diag.size() -1;

    for (int i = 1; i<=num_elements; i++){
        lower_diag[i] = n_mob[i+1]*B_n2[i+1];
    }
}


void Continuity_n::set_rhs(const std::vector<double> &n_mob, const std::vector<double> &B_n1, const std::vector<double> &B_n2, double n_leftBC, double n_rightBC, std::vector<double> &Un){
    for(int i = 1; i<=rhs.size() -1; i++){
        rhs[i] = -Cn*Un[i];
    }
    //BCs
    rhs[1] -= n_mob[0]*B_n2[1]*n_leftBC;
    rhs[rhs.size()-1] -= n_mob[rhs.size()]*B_n1[rhs.size()]*n_rightBC;
}


