#ifndef PARAMETERS_H
#define PARAMETERS_H

//Physical Constants (these should not be changed)
const double q =  1.60217646e-19;         //elementary charge, C
const double kb = 1.3806503e-23;          //Boltzmann const., J/k
const double T = 296.;                      //temperature
const double Vt = (kb*T)/q;                  //thermal voltage
const double epsilon_0 =  8.85418782e-12; //F/m

//Thicknesses (in m)
const double L = 300.0e-9;             //device length in meters. Note: if change this and using a gen_rate file, need to update the number of elements in file to = num_cell-2
const double dx = 1.e-9;
const int num_cell = L/dx;            //number of cells

const double N_LUMO = 10e24;
const double N_HOMO = 10e24;
const double N = 10e24;     //scaling factor helps CV be on order of 1
const double Nsqrd = N*N;

const double G_max =  7e27;  //max photogeneration rate

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
const double E_trap = active_VB + E_gap/2.0;  //this is assumption used--> trap assisted recombo is most effective when trap is located mid-gap--> take VB and add 1/2 of bandgap
const double n1 = N_LUMO*exp(-(active_CB - E_trap)/Vt);
const double p1 = N_HOMO*exp(-(E_trap - active_VB)/Vt);

//////////////////////////////////////////
//Voltage sweep setup
const double Va_min = -0.5;
const double Va_max = 1.2;
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
