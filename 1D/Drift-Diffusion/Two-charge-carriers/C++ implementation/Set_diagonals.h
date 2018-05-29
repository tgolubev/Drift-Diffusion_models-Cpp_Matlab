#ifndef SET_AV_DIAGS_H
#define SET_AV_DIAGS_H

#include<vector>
#include "parameters.h"

void set_main_AVdiag(const std::vector<double> &epsilon, std::vector<double> &a);
void  set_upper_AVdiag(const std::vector<double> &epsilon, std::vector<double> &b);
void  set_lower_AVdiag(const std::vector<double> &epsilon, std::vector<double> &c);
void  set_main_An_diag(const std::vector<double> &n_mob, const std::vector<double> &B_n1, const std::vector<double> &B_n2, std::vector<double> &main_diag);
void  set_upper_An_diag(const std::vector<double> &n_mob, const std::vector<double> &B_n1, std::vector<double> &upper_diag);
void  set_lower_An_diag(const std::vector<double> &n_mob, const std::vector<double> &B_n2, std::vector<double> &lower_diag);
void set_main_Ap_diag(const std::vector<double> &p_mob, const std::vector<double> &B_p1, const std::vector<double> &B_p2, std::vector<double> &main_diag);
void set_upper_Ap_diag(const std::vector<double> &p_mob, const std::vector<double> &B_p2, std::vector<double> &upper_diag);
void set_lower_Ap_diag(const std::vector<double> &p_mob, const std::vector<double> &B_p1, std::vector<double> &lower_diag);

#endif // SET_AV_DIAGS_H

