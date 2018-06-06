#include<vector>
#include"parameters.h"
#include <iostream>

void set_bp(const std::vector<double> &p_mob,const  std::vector<double> &B_p1, const std::vector<double> &B_p2, double p_leftBC, double p_rightBC, std::vector<double> &Up, std::vector<double> &bp){
    for(int i = 1; i<=num_cell-1; i++){
        bp[i] = -Cp*Up[i];
    }
    //BCs
    bp[1] -= p_mob[0]*B_p1[1]*p_leftBC;
    bp[num_cell-1] -= p_mob[num_cell]*B_p2[num_cell]*p_rightBC;
}

void set_bn(const std::vector<double> &n_mob, const std::vector<double> &B_n1, const std::vector<double> &B_n2, double n_leftBC, double n_rightBC, std::vector<double> &Un, std::vector<double> &bn){
    for(int i = 1; i<=num_cell-1; i++){
        bn[i] = -Cn*Un[i];
    }
    //BCs
    bn[1] -= n_mob[0]*B_n2[1]*n_leftBC;
    bn[num_cell-1] -= n_mob[num_cell]*B_n1[num_cell]*n_rightBC;
}

void set_bV(const std::vector<double> &epsilon, const std::vector<double> &n, const std::vector<double> &p, double V_leftBC, double V_rightBC, std::vector<double> &bV){
    for(int i = 1;i<=num_cell-1; i++){
        bV[i] = CV*(n[i] - p[i]);
        //std::cout << "bV " << bV[i] << std::endl;
    }
    bV[1] -= epsilon[0]*V_leftBC;
    bV[num_cell-1] -= epsilon[num_cell]*V_rightBC;
}
