#ifndef CONTINUITY_P_H
#define CONTINUITY_P_H

#include <vector>
#include "simulation.h"  //needs this to know what Simulation is
#include "parameters.h"

class Continuity_p
{
public:
    Continuity_p(Simulation &simul): main_diag(simul.get_num_cell()), upper_diag(simul.get_num_cell()-1), lower_diag(simul.get_num_cell()-1), rhs(simul.get_num_cell()) { }

    void set_main_diag(const std::vector<double> &p_mob,const  std::vector<double> &B_p1,const  std::vector<double> &B_p2);
    void set_upper_diag(const std::vector<double> &p_mob,const  std::vector<double> &B_p1);
    void set_lower_diag(const std::vector<double> &p_mob,const  std::vector<double> &B_p2);
    void set_rhs(const std::vector<double> &p_mob,const  std::vector<double> &B_p1, const std::vector<double> &B_p2, double p_leftBC, double p_rightBC, std::vector<double> &Up);

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
    const double Cp = dx*dx/(Vt*N*mobil);  //can't use static, b/c dx wasn't defined as const, so at each initialization of Continuity_p object, new const will be made.
};

#endif // CONTINUITY_P_H
