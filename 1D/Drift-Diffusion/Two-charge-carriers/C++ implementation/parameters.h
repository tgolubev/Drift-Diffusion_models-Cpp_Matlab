#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>

//Physical Constants (these should not be changed)  (constexpr b/c known at compile-time)
constexpr double q =  1.60217646e-19;         //elementary charge, C
constexpr double kb = 1.3806503e-23;          //Boltzmann const., J/k
constexpr double T = 296.;                      //temperature
constexpr double Vt = (kb*T)/q;                  //thermal voltage
constexpr double epsilon_0 =  8.85418782e-12; //F/m

const double N_LUMO = 10e24;
const double N_HOMO = 10e24;
const double N = N_HOMO;     //scaling factor helps CV be on order of 1
const double Nsqrd = N*N;

const double Photogen_scaling =  7e27;  //max photogeneration rate

//injection barriers
const double phi_a = 0.2;
const double phi_c = 0.1;

//Relative dielectric constants (epsilon/epsilon_0)
const double eps_active = 3.0;

//Mobilities
const double p_mob_active =  4.5e-6;
const double n_mob_active =  4.5e-6;

const double mobil = 5e-6; //scaling for mobility

//Energetics
const double E_gap = 1.5;        //active layer bandgap (eV)
const double active_CB = -3.9;   //active layer conduction band energy (eV)
const double active_VB = -5.4;   //active layer valence band energy (eV)

const double WF_anode = 4.8;
const double WF_cathode = 3.7;

const double Vbi = WF_anode - WF_cathode +phi_a +phi_c;

//Recombination parameters
const double k_rec= 6e-17;  //m^3/s  Langevin recombination for active layer

//////////////////////////////////////////


const double dx = 1.e-9;



#endif // PARAMETERS_H
