#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "constants.h"


struct Parameters   //parameters need to be accessble, so all members are public.
{
    Parameters();
    void reduce_w(){w = w/w_reduce_factor;}
    void relax_tolerance(){tolerance = tolerance*tol_relax_factor;}
    void use_tolerance_eq() {tolerance = tolerance_eq;}
    void use_tolerance_i() {tolerance = tolerance_i;}
    void use_w_i() {w = w_i;}
    void use_w_eq() {w = w_eq;}
    void Initialize();

    //getters
    double get_tolerance() const {return tolerance;}
    double get_w() const {return w;}   //const specifies that fnc shouldn't change w
    int get_num_cell() const {return num_cell;}
    //int get_num_V() const {return num_V;}
    double get_Va_min() const {return Va_min;}
    double get_Va_max() const {return Va_max;}
    double get_increment() const {return increment;}
    //double dx() const {return dx;}
    double N_LUMO, N_HOMO, phi_a, phi_c, eps_active, p_mob_active, n_mob_active;
    double dx, mobil;
    double E_gap, active_CB, active_VB, WF_anode, WF_cathode, N, Nsqrd;
    double Photogen_scaling, k_rec_input;

    //note: seems CANNOT declare unitiliazed constants here and then initialize them later. IF IT IS A CONST, then you must specify the value
    //when it is declared.
    //if want to declare+define/initialize const in a class: use 'static const' or 'static constexpr' (if known at compile time) --> makes sure there is just 1 per class, instead of 1 per object in the class
    double w;
    double w_reduce_factor;
    double tolerance;
    double tol_relax_factor;
    double Vmin, Vmax;
    double tolerance_eq, tolerance_i, w_i, w_eq;
    double L;
    int num_cell;

    double Va_min, Va_max, increment;


};

#endif // PARAMETERS_H
