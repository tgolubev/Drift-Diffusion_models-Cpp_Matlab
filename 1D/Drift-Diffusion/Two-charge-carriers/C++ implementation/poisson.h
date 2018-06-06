#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include "simulation.h"  //needs this to know what Simulation is

class Poisson
{
public:
    Poisson(Simulation &simul): main_diag(simul.get_num_cell()), upper_diag(simul.get_num_cell()-1), lower_diag(simul.get_num_cell()-1), rhs(simul.get_num_cell()) { } //constructor

    void set_main_diag(const std::vector<double> &epsilon);
    void set_upper_diag(const std::vector<double> &epsilon);
    void set_lower_diag(const std::vector<double> &epsilon);
    void set_rhs(const std::vector<double> &epsilon, const std::vector<double> &n, const std::vector<double> &p, double V_leftBC, double V_rightBC);

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
};

#endif // POISSON_H
