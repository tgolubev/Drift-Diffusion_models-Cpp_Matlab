#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include <Eigen/Dense>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include<Eigen/SparseQR>
#include <Eigen/OrderingMethods>

#include "parameters.h"  //needs this to know what paramsation is
#include "constants.h"

class Poisson
{
public:
    Poisson(const Parameters &params);

    //!Setup Poisson equation matrix AV. Note that the matrix will usually be the same throughout the simulation
    //! since it only depends on dielectric constant, so this function should be called only once.
    void setup_matrix();

    //!Setup the right hand side of Poisson equation. This depends on the electron density \param n_matrix,
    //! hole density \param p_matrix, and left and right boundary conditions \param V_leftBC and \param V_rightBC
    void set_rhs(const Eigen::MatrixXd n_matrix, const Eigen::MatrixXd p_matrix);

    void Poisson::to_matrix(const std::vector<double> &V);

    //setters for BC's:
    //for left and right BC's, will use input from the n matrix to determine
    void set_V_topBC(const Parameters &params, double Va);  //need applied voltage input, to know what the BC's are
    void set_V_bottomBC(const Parameters &params, double Va);
    void set_V_leftBC(const std::vector<double> &V);
    void set_V_rightBC(const std::vector<double> &V);

    //getters
    Eigen::VectorXd get_rhs() const {return VecXd_rhs;}  //returns the Eigen object
    Eigen::SparseMatrix<double> get_sp_matrix() const {return sp_matrix;}
    std::vector<double> get_V_topBC() const {return V_topBC;}    //top and bottom  bc getters are needed to determine initial V
    std::vector<double> get_V_bottomBC() const {return V_bottomBC;}
    Eigen::MatrixXd get_V_matrix() const {return V_matrix;}

    //The below getters can be useful for testing and debugging
    //std::vector<double> get_main_diag() const {return main_diag;}
    //std::vector<double> get_upper_diag() const {return upper_diag;}
    //std::vector<double> get_lower_diag() const {return lower_diag;}
    //std::vector<double> get_far_lower_diag() const {return far_lower_diag;}
    //std::vector<double> get_far_upper_diag() const {return far_upper_diag;}
    //std::vector<double> get_V_leftBC() const {return V_leftBC;}
    //std::vector<double> get_V_rightBC() const {return V_rightBC;}

private:
    double CV;    //Note: relative permitivity was moved into the matrix
    int N;  //for convenience define this --> is the number of points in 1D inside the device
    int num_elements;  //for convience so don't have to keep writing params.
    int num_cell;

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
    Eigen::VectorXd VecXd_rhs;  //rhs in Eigen object vector form, for sparse matrix solver
    Eigen::SparseMatrix<double> sp_matrix;
    Eigen::MatrixXd V_matrix;
    Eigen::MatrixXd netcharge;

    //Boundary conditions
    std::vector<double> V_leftBC, V_rightBC, V_bottomBC, V_topBC;

    //!This matrix stores the possibly position-dependent relative dielectric constant.
    Eigen::MatrixXd epsilon;

};

#endif // POISSON_H
