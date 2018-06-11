#include "continuity_n.h"

Continuity_n::Continuity_n(Parameters &params)
{
   main_diag.resize(params.num_cell);
   upper_diag.resize(params.num_cell-1);
   lower_diag.resize(params.num_cell-1);
   rhs.resize(params.num_cell);
   n_mob.resize(params.num_cell+1);

   Cn = params.dx*params.dx/(Vt*params.N*params.mobil);
   n_leftBC = (params.N_LUMO*exp(-(params.E_gap - params.phi_a)/Vt))/params.N;       //this is anode
   n_rightBC = (params.N_LUMO*exp(-params.phi_c/Vt))/params.N;

   std::fill(n_mob.begin(), n_mob.end(), params.n_mob_active/params.mobil);
}



void Continuity_n::setup_eqn(std::vector<double> &B_n1, std::vector<double> &B_n2, std::vector<double> &Un)
{
    set_main_diag(B_n1, B_n2);
    set_upper_diag(B_n1);
    set_lower_diag(B_n2);
    set_rhs(B_n1, B_n2, Un);
}

//-------------------------------Setup An diagonals (Continuity/drift-diffusion solve)-----------------------------
void Continuity_n::set_main_diag(const  std::vector<double> &B_n1,const  std::vector<double> &B_n2)
{
    for (int i = 1; i < main_diag.size(); i++) {
        main_diag[i] = -(n_mob[i]*B_n1[i] + n_mob[i+1]*B_n2[i+1]);
    }
}

//this is b in tridiag_solver
void Continuity_n::set_upper_diag(const  std::vector<double> &B_n1)
{
    for (int i = 1; i < upper_diag.size(); i++) {
        upper_diag[i] = n_mob[i+1]*B_n1[i+1];
    }
}

//this is c in tridiag_solver
void Continuity_n::set_lower_diag(const  std::vector<double> &B_n2)
{
    for (int i = 1; i < lower_diag.size(); i++) {
        lower_diag[i] = n_mob[i+1]*B_n2[i+1];
    }
}


void Continuity_n::set_rhs(const std::vector<double> &B_n1, const std::vector<double> &B_n2, std::vector<double> &Un)
{
    for (int i = 1; i < rhs.size(); i++) {
        rhs[i] = -Cn*Un[i];
    }
    //BCs
    rhs[1] -= n_mob[0]*B_n2[1]*n_leftBC;
    rhs[rhs.size()-1] -= n_mob[rhs.size()]*B_n1[rhs.size()]*n_rightBC;
}


