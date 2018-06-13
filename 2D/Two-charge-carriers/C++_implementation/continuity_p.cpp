#include "continuity_p.h"

Continuity_p::Continuity_p(const Parameters &params)
{
    num_cell = params.num_cell;
    main_diag.resize(num_cell);
    upper_diag.resize(num_cell-1);
    lower_diag.resize(num_cell-1);
    rhs.resize(num_cell);

    Bp_posX.resize(num_cell+1, num_cell+1);  //allocate memory for the matrix object
    Bp_negX.resize(num_cell+1, num_cell+1);
    Bp_posZ.resize(num_cell+1, num_cell+1);
    Bp_negZ.resize(num_cell+1, num_cell+1);

    p_mob.resize(num_cell+1, num_cell+1);  //note: p_mob is an Eigen Matrix object...

    Cp = params.dx*params.dx/(Vt*params.N_dos*params.mobil);  //can't use static, b/c dx wasn't defined as const, so at each initialization of Continuity_p object, new const will be made.

    for (int j =  0; j <= N; j++) {
        p_bottomBC[j] = params.N_HOMO*exp(-params.phi_a/Vt)/params.N_dos;
        p_topBC[j] = params.N_HOMO*exp(-(params.E_gap-params.phi_c)/Vt)/params.N_dos;
    }

    num_elements = params.num_elements;
    N = params.num_cell - 1;
}

//Calculates Bernoulli fnc values, then sets the diagonals and rhs
void Continuity_p::setup_eqn(const Eigen::MatrixXd &V_matrix, const Eigen::MatrixXd &Up_matrix)
{
    Bernoulli_p_X(V_matrix);
    Bernoulli_p_Z(V_matrix);
    set_main_diag();
    set_upper_diag();
    set_lower_diag();
    set_rhs(Up_matrix);
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


void Continuity_p::set_rhs(const Eigen::MatrixXd &Up_matrix)
{
    for (int i = 1; i < rhs.size(); i++) {
        rhs[i] = -Cp*Up[i];
    }
    //BCs
    rhs[1] -= p_mob[0]*B_p1[1]*p_leftBC;
    rhs[rhs.size()-1] -= p_mob[rhs.size()]*B_p2[rhs.size()]*p_rightBC;
}

//---------------------------

void Continuity_p::Bernoulli_p_X(const Eigen::MatrixXd &V_matrix)
{
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(num_cell+1,num_cell+1);

    for (int i = 1; i < num_cell+1; i++)
       for (int j = 1; j < num_cell+1; j++)
           dV(i,j) =  V_matrix(i,j)-V_matrix(i-1,j);

    for (int i = 1; i < V.size(); i++) {
        for (int j = 1; j < num_cell+1; j++) {
            Bp_posX(i,j) = dV(i,j)/(exp(dV(i,j)) - 1.0);
            Bp_negX(i,j) = Bp_posX(i,j)*exp(dV(i,j));
        }
    }
}

void Continuity_p::Bernoulli_p_Z(const Eigen::MatrixXd &V_matrix)
{
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(num_cell+1,num_cell+1);

    for (int i = 1; i < num_cell+1; i++)
       for (int j = 1; j < num_cell+1; j++)
           dV(i,j) =  V_matrix(i,j)-V_matrix(i,j-1);

    for (int i = 1; i < V.size(); i++) {
        for (int j = 1; j < num_cell+1; j++) {
            Bp_posZ(i,j) = dV(i,j)/(exp(dV(i,j)) - 1.0);
            Bp_negZ(i,j) =  Bp_posZ(i,j)*exp(dV(i,j));
        }
    }
}
