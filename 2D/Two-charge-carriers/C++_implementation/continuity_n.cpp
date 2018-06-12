#include "continuity_n.h"

Continuity_n::Continuity_n(const Parameters &params)
{
   num_cell = params.num_cell;
   main_diag.resize(num_cell);
   upper_diag.resize(num_cell-1);
   lower_diag.resize(num_cell-1);
   rhs.resize(num_cell);
   B_n1.resize(num_cell+1);
   B_n2.resize(num_cell+1);
   //n_mob.resize(params.num_cell +1, std::vector<double>(params.num_cell+1));

   n_mob.resize(num_cell+1, num_cell+1);  //note: n_mob is an Eigen Matrix object...
/*
 //for now  fill matrix, by using fill for each of the vectors inside

   for (int i = 0; i < params.num_cell; i++) {
       std::fill(n_mob[i].begin(), n_mob[i].end(), params.n_mob_active/params.mobil);
   }
   */

   Cn = params.dx*params.dx/(Vt*params.N_dos*params.mobil);

   n_bottomBC = params.N_LUMO*exp(-(params.E_gap-params.phi_a)/Vt)/params.N_dos;
   n_topBC = params.N_LUMO*exp(-params.phi_c/Vt)/params.N_dos;



}

//Calculates Bernoulli fnc values, then sets the diagonals and rhs
//use the V_matrix for setup, to be able to write equations in terms of (x,z) coordingates
void Continuity_n::setup_eqn(const Eigen::MatrixXd &V_matrix, const Eigen::MatrixXd &Un_matrix)
{
    BernoulliFnc_n(V_matrix);
    set_main_diag();
    set_upper_diag();
    set_lower_diag();
    set_rhs(Un_matrix);

}

//-------------------------------Setup An diagonals (Continuity/drift-diffusion solve)-----------------------------
void Continuity_n::set_main_diag()
{
    for (int i = 1; i < main_diag.size(); i++) {
        main_diag[i] = -(n_mob[i]*B_n1[i] + n_mob[i+1]*B_n2[i+1]);
    }
}

//this is b in tridiag_solver
void Continuity_n::set_upper_diag()
{
    for (int i = 1; i < upper_diag.size(); i++) {
        upper_diag[i] = n_mob[i+1]*B_n1[i+1];
    }
}

//this is c in tridiag_solver
void Continuity_n::set_lower_diag()
{
    for (int i = 1; i < lower_diag.size(); i++) {
        lower_diag[i] = n_mob[i+1]*B_n2[i+1];
    }
}


void Continuity_n::set_rhs(const Eigen::MatrixXd &Un_matrix)
{
    for (int i = 1; i < rhs.size(); i++) {
        rhs[i] = -Cn*Un[i];
    }
    //BCs
    rhs[1] -= n_mob[0]*B_n2[1]*n_leftBC;
    rhs[rhs.size()-1] -= n_mob[rhs.size()]*B_n1[rhs.size()]*n_rightBC;
}

//------------------------
//Note: are using the V matrix for Bernoulli calculations.
//Makes it clearer to write indices in terms of (x,z) real coordinate values.

void Continuity_n::Bernoulli_n_X(const Eigen::MatrixXd &V_matrix)
{
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(num_cell+1,num_cell+1);

    for (int i = 1; i < num_cell+1; i++)             //Note: the indices are shifted by 1 from Matlab version (here bndry is at index 0)
        for (int j = 1; j < num_cell+1; j++)
            dV(i,j) =  V_matrix(i,j)-V_matrix(i-1,j);


    for (int i = 1; i < num_cell+1; i++) {           //note: the indexing done a bit different than Matlab (see 1D C++ version)
        for (int j = 1; j < num_cell+1; j++) {
            Bn_posX(i,j) = dV(i,j)/(exp(dV(i,j)) - 1.0);
            Bn_negX(i,j) = Bn_posX(i,j)*exp(dV(i,j));
        }
    }
}

void Continuity_n::Bernoulli_n_Z(const Eigen::MatrixXd &V_matrix)
{
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(num_cell+1,num_cell+1);

    for (int i = 1; i < num_cell+1; i++)
       for (int j = 1; j < num_cell+1; j++)
           dV(i,j) =  V_matrix(i,j)-V_matrix(i,j-1);

    for (int i = 1; i < num_cell+1; i++) {
        for (int j = 1; j < num_cell+1; j++) {
            Bn_posZ(i,j) = dV(i,j)/(exp(dV(i,j)) - 1.0);
            Bn_negZ(i,j) =  Bn_posZ(i,j)*exp(dV(i,j));
        }
    }
}
