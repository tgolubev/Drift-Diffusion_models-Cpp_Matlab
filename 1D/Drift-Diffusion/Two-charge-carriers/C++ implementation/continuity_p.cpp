#include "continuity_p.h"

//------------------------------Setup Ap diagonals----------------------------------------------------------------
void Continuity_p::set_main_diag(const std::vector<double> &p_mob,const  std::vector<double> &B_p1, const std::vector<double> &B_p2){
    int num_elements = main_diag.size() - 1;

    for (int i=1; i<= num_elements;i++){
        main_diag[i] = -(p_mob[i]*B_p2[i] + p_mob[i+1]*B_p1[i+1]);
    }
}

//this is b in tridiag_solver
void Continuity_p::set_upper_diag(const std::vector<double> &p_mob,const  std::vector<double> &B_p2){
    int num_elements = upper_diag.size() -1;  //-1 since vectors fill from 0, but I'm indexing from 1
    //NOTE: num_elements here, is actually # of elements in the off-diagonal

    for (int i = 1; i<=num_elements; i++){
        upper_diag[i] = p_mob[i+1]*B_p2[i+1];
    }
}

//this is c in tridiag_solver
void Continuity_p::set_lower_diag(const std::vector<double> &p_mob,const  std::vector<double> &B_p1){
    int num_elements = lower_diag.size() -1;

    for (int i = 1; i<=num_elements; i++){
        lower_diag[i] = p_mob[i+1]*B_p1[i+1];
    }
}


void Continuity_p::set_rhs(const std::vector<double> &p_mob,const  std::vector<double> &B_p1, const std::vector<double> &B_p2, double p_leftBC, double p_rightBC, std::vector<double> &Up){
    for(int i = 1; i<=rhs.size()-1; i++){
        rhs[i] = -Cp*Up[i];
    }
    //BCs
    rhs[1] -= p_mob[0]*B_p1[1]*p_leftBC;
    rhs[rhs.size()-1] -= p_mob[rhs.size()]*B_p2[rhs.size()]*p_rightBC;
}
