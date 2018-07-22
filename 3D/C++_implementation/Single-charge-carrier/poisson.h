#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include <Eigen/Dense>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include<Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include <unsupported/Eigen/CXX11/Tensor>

#include "parameters.h"  //needs this to know what paramsation is
#include "constants.h"

class Poisson
{

typedef Eigen::Triplet<double> Trp;

public:
    Poisson(const Parameters &params);

    //!Setup Poisson equation matrix AV. Note that the matrix will usually be the same throughout the simulation
    //! since it only depends on dielectric constant, so this function should be called only once.
    void setup_matrix();

    //!Setup the right hand side of Poisson equation. This depends on the hole density \param p.
    void set_rhs(const Eigen::Tensor<double, 3> &p);

    //setters for BC's:
    //for left and right BC's, will use input from the n matrix to determine
    void set_V_topBC(const Parameters &params, double Va);  //need applied voltage input, to know what the BC's are
    void set_V_bottomBC(const Parameters &params, double Va);

    //getters
    Eigen::VectorXd get_rhs() const {return VecXd_rhs;}  //returns the Eigen object
    Eigen::SparseMatrix<double> get_sp_matrix() const {return sp_matrix;}
    double get_V_topBC(int i, int j) const {return V_topBC(i,j);}    //top and bottom  bc getters are needed to determine initial V
    double get_V_bottomBC(int i, int j) const {return V_bottomBC(i,j);}
    Eigen::Tensor<double, 3> get_V_matrix() const {return V_matrix;}

    //The below getters can be useful for testing and debugging
    //std::vector<double> get_main_diag() const {return main_diag;}
    //std::vector<double> get_upper_diag() const {return upper_diag;}
    //std::vector<double> get_lower_diag() const {return lower_diag;}
    //std::vector<double> get_far_lower_diag() const {return far_lower_diag;}
    //std::vector<double> get_far_upper_diag() const {return far_upper_diag;}


private:
    double CV;    //Note: relative permitivity was moved into the matrix
    int Nx, Ny, Nz;  //for convenience define this --> is the number of points in 1D inside the device
    int num_elements;  //for convience so don't have to keep writing params.
    int num_cell;

    void set_far_lower_diag();
    void set_lower_diag();
    void set_main_lower_diag();
    void set_main_diag();
    void set_main_upper_diag();
    void set_upper_diag();
    void set_far_upper_diag();

    std::vector<double> far_lower_diag;
    std::vector<double> far_upper_diag;
    std::vector<double> main_lower_diag;
    std::vector<double> main_diag;
    std::vector<double> main_upper_diag;
    std::vector<double> upper_diag;
    std::vector<double> lower_diag;
    std::vector<double> rhs;
    Eigen::VectorXd VecXd_rhs;  //rhs in Eigen object vector form, for sparse matrix solver
    Eigen::SparseMatrix<double> sp_matrix;
    Eigen::Tensor<double, 3> V_matrix;
    Eigen::Tensor<double, 3> netcharge;

    std::vector<Trp> triplet_list;
    int trp_cnt;  //for counting the triplets

    //Boundary conditions
    Eigen::MatrixXd V_bottomBC, V_topBC;

    //!This matrix stores the possibly position-dependent relative dielectric constant.
    Eigen::Tensor<double, 3> epsilon, epsilon_avg_X, epsilon_avg_Y, epsilon_avg_Z;

};

#endif // POISSON_H
