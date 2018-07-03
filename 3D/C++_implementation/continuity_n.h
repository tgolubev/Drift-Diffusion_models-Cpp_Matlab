#ifndef CONTINUITY_N_H
#define CONTINUITY_N_H

#include <vector>
#include <Eigen/Dense>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include<Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include <unsupported/Eigen/CXX11/Tensor>

#include "parameters.h"  //needs this to know what parameters are
#include "constants.h"

class Continuity_n
{

typedef Eigen::Triplet<double> Trp;

public:
    Continuity_n(const Parameters &params);

    //!Sets up the matrix equation An*n = bn for continuity equation for electrons.
    //!\param V stores the voltage and is needed to calculate Bernoulli fnc.'s.
    //!\param Un stores the net generation rate, needed for the right hand side.
    //!\param n the electron density is needed to setup the boundary conditions.
    void setup_eqn(const Eigen::Tensor<double, 3> &V_matrix, const std::vector<double> &Un, const std::vector<double> &n);

    void calculate_currents();

   void Continuity_n::to_matrix(const std::vector<double> &n);

    //setters for BC's:
    //for left and right BC's, will use input from the n matrix to determine
    void set_n_topBC();
    void set_n_bottomBC();
    void set_n_leftBC_X(const std::vector<double> &n);
    void set_n_rightBC_X(const std::vector<double> &n);
    void set_n_leftBC_Y(const std::vector<double> &n);
    void set_n_rightBC_Y(const std::vector<double> &n);

    //getters (const keyword ensures that fnc doesn't change anything)
    Eigen::VectorXd get_rhs() const {return VecXd_rhs;}  //returns the Eigen object
    Eigen::SparseMatrix<double> get_sp_matrix() {return sp_matrix;}
    Eigen::Tensor<double, 3> get_n_matrix() const {return n_matrix;}
    Eigen::MatrixXd get_n_bottomBC() const {return n_bottomBC;}  //bottom and top are needed to set initial conditions
    Eigen::MatrixXd get_n_topBC() const {return n_topBC;}

    Eigen::Tensor<double, 3> get_Jn_X() const {return Jn_X;}
    Eigen::Tensor<double, 3> get_Jn_Y() const {return Jn_Y;}
    Eigen::Tensor<double, 3> get_Jn_Z() const {return Jn_Z;}

    //The below getters can be useful for testing and debugging
    //std::vector<double> get_main_diag() const {return main_diag;}
    //std::vector<double> get_upper_diag() const {return upper_diag;}
    //std::vector<double> get_lower_diag() const {return lower_diag;}
    //std::vector<double> get_far_upper_diag() const {return far_upper_diag;}
    //std::vector<double> get_far_lower_diag() const {return far_lower_diag;}
    //Eigen::MatrixXd  get_Bn_posX() const {return Bn_posX;}
    //Eigen::MatrixXd  get_Bn_negX() const {return Bn_negX;}
    //Eigen::MatrixXd  get_Bn_posZ() const {return Bn_posZ;}
    //Eigen::MatrixXd  get_Bn_negZ() const {return Bn_negZ;}
    //std::vector<double> get_n_leftBC() const {return n_leftBC;}
    //std::vector<double> get_n_rightBC() const {return n_rightBC;}
    //Eigen::MatrixXd get_n_mob() const {return n_mob;}

private:
    std::vector<double> far_lower_diag;
    std::vector<double> far_upper_diag;
    std::vector<double> main_lower_diag;
    std::vector<double> main_diag;
    std::vector<double> main_upper_diag;
    std::vector<double> upper_diag;
    std::vector<double> lower_diag;
    std::vector<double> rhs;
    Eigen::Tensor<double, 3> n_mob;  //!Matrix storing the position dependent electron mobility
    Eigen::Tensor<double, 3> n_mob_avg_X, n_mob_avg_Y, n_mob_avg_Z;
    Eigen::VectorXd VecXd_rhs;  //rhs in Eigen object vector form, for sparse matrix solver
    Eigen::SparseMatrix<double> sp_matrix;
    Eigen::Tensor<double, 3> n_matrix;
    Eigen::Tensor<double, 3> Jn_Z;
    Eigen::Tensor<double, 3> Jn_X;
    Eigen::Tensor<double, 3> Jn_Y;

    Eigen::Tensor<double, 3> values;  //temporary variable to store values for diagonal filling
    Eigen::Tensor<double, 3> values2, values3, values4;  //these additional temp variables needed for main diag filling

    std::vector<Trp> triplet_list;
    int trp_cnt;  //for counting the triplets
    double J_coeff;  //coefficient for curents eqn

    //Boundary conditions
    Eigen::MatrixXd n_leftBC_X, n_rightBC_X, n_leftBC_Y, n_rightBC_Y, n_bottomBC, n_topBC;

    Eigen::Tensor<double, 3> Bn_posX;  //bernoulli (+dV_x)
    Eigen::Tensor<double, 3> Bn_negX;  //bernoulli (-dV_x)
    Eigen::Tensor<double, 3> Bn_posY;  //bernoulli (+dV_y)
    Eigen::Tensor<double, 3> Bn_negY;  //bernoulli (-dV_y)
    Eigen::Tensor<double, 3> Bn_posZ;  //bernoulli (+dV_z)
    Eigen::Tensor<double, 3> Bn_negZ;  //bernoulli (-dV_z)

    double Cn;
    int num_cell, num_elements;
    int N;

    //!Calculates the Bernoulli functions for dV in x direction and updates member arrays
    void Bernoulli_n_X(const Eigen::Tensor<double, 3> &V_matrix);

    //!Calculates the Bernoulli functions for dV in y direction and updates member arrays
    void Bernoulli_n_Y(const Eigen::Tensor<double, 3> &V_matrix);

    //!Calculates the Bernoulli functions for dV in z direction and updates member arrays
    void Bernoulli_n_Z(const Eigen::Tensor<double, 3> &V_matrix);

    //matrix setup functions
    void set_far_lower_diag();
    void set_lower_diag();
    void set_main_lower_diag();
    void set_main_diag();
    void set_main_upper_diag();
    void set_upper_diag();
    void set_far_upper_diag();
    void set_rhs(const std::vector<double> &Un);
};

#endif // CONTINUITY_N_H
