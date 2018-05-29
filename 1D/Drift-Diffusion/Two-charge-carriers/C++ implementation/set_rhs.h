#ifndef SET_RHS_H
#define SET_RHS_H

#include<vector>
#include "parameters.h"

void set_bp(const std::vector<double> &p_mob,const  std::vector<double> &B_p1, const std::vector<double> &B_p2, double p_leftBC, double p_rightBC, std::vector<double> &Up, std::vector<double> &bp);
void set_bn(const std::vector<double> &n_mob,const  std::vector<double> &B_n1, const std::vector<double> &B_n2, double n_leftBC, double n_rightBC, std::vector<double> &Un, std::vector<double> &bn);
void set_bV(const std::vector<double> &epsilon,const  std::vector<double> &n, const std::vector<double> &p, double V_leftBC, double V_rightBC, std::vector<double> &bV);

#endif // SET_RHS_H
