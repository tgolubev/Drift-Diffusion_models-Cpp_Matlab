#ifndef CONTINUITY_N_H
#define CONTINUITY_N_H

#include <vector>
#include "simulation.h"  //needs this to know what Simulation is
#include "parameters.h"

class Continuity_n
{
public:
    Continuity_n(Simulation &simul): main_diag(simul.get_num_cell()), upper_diag(simul.get_num_cell()-1), lower_diag(simul.get_num_cell()-1),  rhs(simul.get_num_cell()) { }

    void set_main_diag(const std::vector<double> &n_mob,const  std::vector<double> &B_n1,const  std::vector<double> &B_n2);
    void set_upper_diag(const std::vector<double> &n_mob,const  std::vector<double> &B_n1);
    void set_lower_diag(const std::vector<double> &n_mob,const  std::vector<double> &B_n2);
    void set_rhs(const std::vector<double> &n_mob, const std::vector<double> &B_n1, const std::vector<double> &B_n2, double n_leftBC, double n_rightBC, std::vector<double> &Un);

    //getters
    std::vector<double> get_main_diag() const {return main_diag;}
    std::vector<double> get_upper_diag() const {return upper_diag;}
    std::vector<double> get_lower_diag() const {return lower_diag;}
    std::vector<double> get_rhs() const {return rhs;}

private:
    std::vector<double> main_diag;
    std::vector<double> upper_diag;
    std::vector<double> lower_diag;
    std::vector<double> rhs;
    const double Cn = dx*dx/(Vt*N*mobil);
};

#endif // CONTINUITY_N_H
