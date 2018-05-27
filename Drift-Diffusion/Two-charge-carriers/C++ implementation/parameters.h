#ifndef PARAMETERS_H
#define PARAMETERS_H

//Physical Constants (these should not be changed)
//for now, define all parameters here, in main file, and use extern in parameters.h with declarations--> so then all my files can see all the parameters
//Physical Constants (these should not be changed)
const double q =  1.60217646e-19;         //elementary charge, C
const double kb = 1.3806503e-23;          //Boltzmann const., J/k
const double T = 296.;                      //temperature
const double Vt = (kb*T)/q;                  //thermal voltage
const double epsilon_0 =  8.85418782e-12; //F/m

//Thicknesses (in m)
const double L_HTL = 20e-9;
const double L_active = 100e-9;
const double L_ETL = 20.0e-9;

const double L = 100.0e-9;             //device length in meters:
const double dx = 1.e-9;
const int num_cell = L/dx;            //number of cells   --

const int l_HTL_int = L_HTL/dx;               // position of HTL/active layer interface
const int l_ETL_int = (L_HTL + L_active)/dx;         //position of activeskite/ETL interface

const double N_LUMO = 8.1e24;
const double N_HOMO = 8.1e24;
const double N = 8.1e24;     //scaling factor helps CV be on order of 1
const double Nsqrd = N*N;

const double G_max =  3e27;  //max photogeneration rate


//injection barriers
const double phi_a = 0.2;
const double phi_c = 0.2;

//Dielectric constants
const double eps_active = 3.0;

//Mobilities
const double p_mob_active =  4.5e-6;
const double n_mob_active =  4.5e-6;

const double mobil = 5e-6; //scaling for mobility

//Energetics
const double E_gap = 1.5;        //active layer bandgap

const double WF_anode = 4.8;
const double WF_cathode = 4.2;

const double Vbi = WF_anode - WF_cathode +phi_a +phi_c;

//Recombination parameters
const double k_rec= 6e-17;  //m^3/s  Langevin recombination for activeskite

//////////////////////////////////////////
//Voltage sweep setup
const double Va_min = -0.5;
const double Va_max = 1.3;
const double increment = 0.01;
const int num_V = floor((Va_max-Va_min)/increment)+1;

//Simulation parameters
const double w_eq = 0.01;                     //for 1st equil convergence weighting factor
const double w_i = 0.2;                       // %set up of weighting factor
const double tolerance_eq = 100*5e-12;        //error tolerance for equil run
const double tolerance_i = 5e-12;             //initial error tolerance for 1st non-equil Va (1s tVa with generation)d\

const double Cp = dx*dx/(Vt*N*mobil);          //note: scaled p_mob and n_mob are inside matrices
const double Cn = dx*dx/(Vt*N*mobil);
const double CV = N*dx*dx*q/(epsilon_0*Vt);    //relative permitivity was moved into the matrix

#endif // PARAMETERS_H
