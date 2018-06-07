
#include "bernoulli.h"


void BernoulliFnc_n(const std::vector<double> &V, std::vector<double> &B_n1, std::vector<double> &B_n2)
{

    std::vector<double> dV(V.size());

    for(int i = 1; i<V.size();i++){
        dV[i] =  V[i]-V[i-1];
    }

    for(int i = 1; i<V.size();i++){
        B_n1[i] = dV[i]/(exp(dV[i])-1.0);
        B_n2[i] = B_n1[i]*exp(dV[i]);
    }
}


void BernoulliFnc_p(const std::vector<double> &V, std::vector<double> &B_p1, std::vector<double> &B_p2)
{
    std::vector<double> dV(V.size());

    for(int i = 1; i<V.size();i++){
        dV[i] =  V[i]-V[i-1];
    }

    for(int i = 1; i<V.size();i++){
        B_p1[i] = dV[i]/(exp(dV[i])-1.0);
        B_p2[i] = B_p1[i]*exp(dV[i]);
    }
}
