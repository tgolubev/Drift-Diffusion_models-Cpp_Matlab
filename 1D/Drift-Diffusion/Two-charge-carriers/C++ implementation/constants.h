#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

//Physical Constants (these should not be changed)  (constexpr b/c known at compile-time)
constexpr double q =  1.60217646e-19;         //elementary charge, C
constexpr double kb = 1.3806503e-23;          //Boltzmann const., J/k
constexpr double T = 296.;                    //temperature K
constexpr double Vt = (kb*T)/q;               //thermal voltage V
constexpr double epsilon_0 =  8.85418782e-12; //F/m


#endif // CONSTANTS_H
