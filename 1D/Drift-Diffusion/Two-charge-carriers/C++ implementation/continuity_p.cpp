#include "continuity_p.h"


void Continuity_p::setup_eqn(std::vector<double> &p_mob, std::vector<double> &B_p1, std::vector<double> & B_p2, std::vector<double> & Up)
{
    set_main_diag(p_mob, B_p1, B_p2);
    set_upper_diag(p_mob, B_p2);
    set_lower_diag(p_mob, B_p1);
    set_rhs(p_mob, B_p1, B_p2, Up);
}

//------------------------------Setup Ap diagonals----------------------------------------------------------------
void Continuity_p::set_main_diag(const std::vector<double> &p_mob,const  std::vector<double> &B_p1, const std::vector<double> &B_p2)
{
    for (int i=1; i<main_diag.size();i++){
        main_diag[i] = -(p_mob[i]*B_p2[i] + p_mob[i+1]*B_p1[i+1]);
    }
}

//this is b in tridiag_solver
void Continuity_p::set_upper_diag(const std::vector<double> &p_mob,const  std::vector<double> &B_p2)
{
    for (int i = 1; i<upper_diag.size(); i++){
        upper_diag[i] = p_mob[i+1]*B_p2[i+1];
    }
}

//this is c in tridiag_solver
void Continuity_p::set_lower_diag(const std::vector<double> &p_mob,const  std::vector<double> &B_p1){

    for (int i = 1; i<lower_diag.size(); i++){
        lower_diag[i] = p_mob[i+1]*B_p1[i+1];
    }
}


void Continuity_p::set_rhs(const std::vector<double> &p_mob,const  std::vector<double> &B_p1, const std::vector<double> &B_p2, std::vector<double> &Up)
{
    for(int i = 1; i<rhs.size(); i++){
        rhs[i] = -Cp*Up[i];
    }
    //BCs
    rhs[1] -= p_mob[0]*B_p1[1]*p_leftBC;
    rhs[rhs.size()-1] -= p_mob[rhs.size()]*B_p2[rhs.size()]*p_rightBC;
}
