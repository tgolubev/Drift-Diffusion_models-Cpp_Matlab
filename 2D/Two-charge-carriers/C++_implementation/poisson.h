#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include <Eigen/Dense>

#include "parameters.h"  //needs this to know what paramsation is
#include "constants.h"

class Poisson
{
public:
    Poisson(const Parameters &params);

    //!Setup Poisson equation matrix AV. Note that the matrix will usually be the same throughout the simulation
    //! since it only depends on dielectric constant, so this function should be called only once.
    void setup_matrix();

    //!Setup the right hand side of Poisson equation. This depends on the electron density \param n,
    //! hole density \param p, and left and right boundary conditions \param V_leftBC and \param V_rightBC
    void set_rhs(const std::vector<double> &n, const std::vector<double> &p, std::vector<double> &V_leftBC, std::vector<double> &V_rightBC, double V_bottomBC, double V_topBC);

    //getters
    std::vector<double> get_main_diag() const {return main_diag;}
    std::vector<double> get_upper_diag() const {return upper_diag;}
    std::vector<double> get_lower_diag() const {return lower_diag;}
    std::vector<double> get_far_lower_diag() const {return far_lower_diag;}
    std::vector<double> get_far_upper_diag() const {return far_upper_diag;}
    std::vector<double> get_rhs() const {return rhs;}

private:
    void set_far_lower_diag();
    void set_lower_diag();
    void set_main_diag();
    void set_upper_diag();
    void set_far_upper_diag();

    std::vector<double> far_lower_diag;
    std::vector<double> lower_diag;
    std::vector<double> main_diag;
    std::vector<double> upper_diag;
    std::vector<double> far_upper_diag;
    std::vector<double> rhs;

    Eigen::MatrixXd netcharge;  //will store net charge

    //!This matrix stores the possibly position-dependent relative dielectric constant.
    Eigen::MatrixXd epsilon;
    //std::vector<std::vector<double> > epsilon;

    double CV;    //Note: relative permitivity was moved into the matrix

    int N;  //for convenience define this --> is the number of points in 1D inside the device
    int num_elements;  //for convience so don't have to keep writing params.
    int num_cell;

};

#endif // POISSON_H
