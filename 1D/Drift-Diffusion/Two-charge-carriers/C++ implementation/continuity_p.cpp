#include "continuity_p.h"


Continuity_p::Continuity_p(Parameters &params)
{
    main_diag.resize(params.num_cell);
    upper_diag.resize(params.num_cell-1);
    lower_diag.resize(params.num_cell-1);
    rhs.resize(params.num_cell);
    p_mob.resize(params.num_cell+1);

    std::fill(p_mob.begin(), p_mob.end(), params.p_mob_active/params.mobil);

    Cp = params.dx*params.dx/(Vt*params.N*params.mobil);  //can't use static, b/c dx wasn't defined as const, so at each initialization of Continuity_p object, new const will be made.
    p_leftBC = (params.N_HOMO*exp(-params.phi_a/Vt))/params.N;
    p_rightBC = (params.N_HOMO*exp(-(params.E_gap - params.phi_c)/Vt))/params.N;
}


void Continuity_p::setup_eqn(std::vector<double> &B_p1, std::vector<double> &B_p2, std::vector<double> &Up)
{
    set_main_diag(B_p1, B_p2);
    set_upper_diag(B_p2);
    set_lower_diag(B_p1);
    set_rhs(B_p1, B_p2, Up);
}

//------------------------------Setup Ap diagonals----------------------------------------------------------------
void Continuity_p::set_main_diag(const  std::vector<double> &B_p1, const std::vector<double> &B_p2)
{
    for (int i = 1; i < main_diag.size(); i++) {
        main_diag[i] = -(p_mob[i]*B_p2[i] + p_mob[i+1]*B_p1[i+1]);
    }
}

//this is b in tridiag_solver
void Continuity_p::set_upper_diag(const  std::vector<double> &B_p2)
{
    for (int i = 1; i < upper_diag.size(); i++) {
        upper_diag[i] = p_mob[i+1]*B_p2[i+1];
    }
}

//this is c in tridiag_solver
void Continuity_p::set_lower_diag(const  std::vector<double> &B_p1)
{

    for (int i = 1; i < lower_diag.size(); i++) {
        lower_diag[i] = p_mob[i+1]*B_p1[i+1];
    }
}


void Continuity_p::set_rhs(const  std::vector<double> &B_p1, const std::vector<double> &B_p2, std::vector<double> &Up)
{
    for (int i = 1; i < rhs.size(); i++) {
        rhs[i] = -Cp*Up[i];
    }
    //BCs
    rhs[1] -= p_mob[0]*B_p1[1]*p_leftBC;
    rhs[rhs.size()-1] -= p_mob[rhs.size()]*B_p2[rhs.size()]*p_rightBC;
}
