#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "constants.h"
#include <string>
#include<cmath>
#include <fstream>
#include <iostream>
#include <iomanip>


struct Parameters   //parameters need to be accessble, so all members are public.
{
    Parameters(){}
    void reduce_w(){w = w/w_reduce_factor;}
    void relax_tolerance(){tolerance = tolerance*tol_relax_factor;}
    void use_tolerance_eq() {tolerance = tolerance_eq;}
    void use_tolerance_i() {tolerance = tolerance_i;}
    void use_w_i() {w = w_i;}
    void use_w_eq() {w = w_eq;}

    void Initialize();
    void isPositive(double input, const std::string &comment);
    void isPositive(int input, const std::string &comment);
    void isNegative(double input, const std::string &comment);
    void isNegative(int input, const std::string &comment);

    double N_LUMO, N_HOMO, phi_a, phi_c, eps_active, p_mob_active, n_mob_active;
    double mobil;
    double E_gap, active_CB, active_VB, WF_anode, WF_cathode, N_dos, Nsqrd;
    double Photogen_scaling, k_rec;

    double w;
    double w_reduce_factor;
    double tolerance, tolerance_eq;
    double tol_relax_factor;
    double Vmin, Vmax;

    double tolerance_i, w_i, w_eq;
    double Lx, Ly, Lz;
    double dx, dy, dz;
    int Nx, Ny, Nz;
    int num_cell_x, num_cell_y, num_cell_z, num_elements;  //num_elements = (num_cell-1)^3
    std::string GenRateFileName;
    double Va_min, Va_max, increment;
    double Vbi;


};

#endif // PARAMETERS_H
