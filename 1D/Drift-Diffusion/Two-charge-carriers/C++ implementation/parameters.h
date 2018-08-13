#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "constants.h"
#include <string>
#include<cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>


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
    double dx, mobil;
    double E_gap, active_CB, active_VB, WF_anode, WF_cathode, N, Nsqrd;
    double Photogen_scaling, k_rec;

    double w;
    double w_reduce_factor;
    double tolerance, tolerance_eq;
    double tol_relax_factor;
    double Vmin, Vmax;

    double tolerance_i, w_i, w_eq;
    double L;
    int num_cell;
    std::string GenRateFileName;
    double Va_min, Va_max, increment;

    //optimization (auto fitting) parameters
    bool auto_fit;
    int optim_method;
    double fit_tolerance;
    double optim_max_iter;
    std::string exp_data_file_name;
    double Photogen_scaling_min, Photogen_scaling_max;
    double n_mob_active_max, n_mob_active_min;

    bool PSO_Clerc_Kennedy;

    std::vector<double*> vars;  //will store pointers to the parameters/variables which are optimizing
    std::vector<double> vars_min;  //stores the min of the ranges
    std::vector<double> vars_max;
};

#endif // PARAMETERS_H
