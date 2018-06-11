#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include "parameters.h"  //needs this to know what paramsation is
#include "constants.h"

class Poisson
{
public:
    Poisson(Parameters &params): main_diag(params.num_cell), upper_diag(params.num_cell-1), lower_diag(params.num_cell-1), rhs(params.num_cell)
    {
        CV = params.N*params.dx*params.dx*q/(epsilon_0*Vt);
    } //constructor

    void setup_matrix(std::vector<double> &epsilon);
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
    double CV;    //relative permitivity was moved into the matrix

    void set_main_diag(const std::vector<double> &epsilon);
    void set_upper_diag(const std::vector<double> &epsilon);
    void set_lower_diag(const std::vector<double> &epsilon);
};

#endif // POISSON_H
