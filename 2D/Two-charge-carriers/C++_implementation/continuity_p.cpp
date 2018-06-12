#include "continuity_p.h"

Continuity_p::Continuity_p(const Parameters &params)
{
    num_cell = params.num_cell;
    main_diag.resize(num_cell);
    upper_diag.resize(num_cell-1);
    lower_diag.resize(num_cell-1);
    rhs.resize(num_cell);
    B_p1.resize(num_cell+1);
    B_p2.resize(num_cell+1);

    p_mob.resize(num_cell+1, num_cell+1);  //note: p_mob is an Eigen Matrix object...

    //p_mob.resize(params.num_cell +1, std::vector<double>(params.num_cell+1));

    /*
    //for now fill matrix, by using fill for each of the vectors inside
    for (int i = 0; i < params.num_cell; i++) {
        std::fill(p_mob[i].begin(), p_mob[i].end(), params.p_mob_active/params.mobil);
    }
    */

    Cp = params.dx*params.dx/(Vt*params.N*params.mobil);  //can't use static, b/c dx wasn't defined as const, so at each initialization of Continuity_p object, new const will be made.

    p_bottomBC = params.N_HOMO*exp(-params.phi_a/Vt)/params.N_dos;
    p_topBC = params.N_HOMO*exp(-(params.E_gap-params.phi_c)/Vt)/params.N_dos;


}

//Calculates Bernoulli fnc values, then sets the diagonals and rhs
void Continuity_p::setup_eqn(const std::vector<double> &V, const std::vector<double> &Up)
{
    BernoulliFnc_p(V);
    set_main_diag();
    set_upper_diag();
    set_lower_diag();
    set_rhs(Up);
}

//------------------------------Setup Ap diagonals----------------------------------------------------------------
void Continuity_p::set_main_diag()
{
    for (int i = 1; i < main_diag.size(); i++) {
        main_diag[i] = -(p_mob[i]*B_p2[i] + p_mob[i+1]*B_p1[i+1]);
    }
}

//this is b in tridiag_solver
void Continuity_p::set_upper_diag()
{
    for (int i = 1; i < upper_diag.size(); i++) {
        upper_diag[i] = p_mob[i+1]*B_p2[i+1];
    }
}

//this is c in tridiag_solver
void Continuity_p::set_lower_diag()
{
    for (int i = 1; i < lower_diag.size(); i++) {
        lower_diag[i] = p_mob[i+1]*B_p1[i+1];
    }
}


void Continuity_p::set_rhs(const std::vector<double> &Up)
{
    for (int i = 1; i < rhs.size(); i++) {
        rhs[i] = -Cp*Up[i];
    }
    //BCs
    rhs[1] -= p_mob[0]*B_p1[1]*p_leftBC;
    rhs[rhs.size()-1] -= p_mob[rhs.size()]*B_p2[rhs.size()]*p_rightBC;
}

//---------------------------

void Continuity_p::BernoulliFnc_p(const std::vector<double> &V)
{
    std::vector<double> dV(V.size());

    for (int i = 1; i < V.size(); i++) {
        dV[i] =  V[i]-V[i-1];
    }

    for (int i = 1; i < V.size(); i++) {
        B_p1[i] = dV[i]/(exp(dV[i]) - 1.0);
        B_p2[i] = B_p1[i]*exp(dV[i]);
    }
}
