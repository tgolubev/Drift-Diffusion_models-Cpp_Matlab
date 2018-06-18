#ifndef CONTINUITY_N_H
#define CONTINUITY_N_H

#include <vector>
#include <Eigen/Dense>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include<Eigen/SparseQR>
#include <Eigen/OrderingMethods>

#include "parameters.h"  //needs this to know what parameters are
#include "constants.h"

class Continuity_n
{

//typedef Eigen::Triplet<double> Trp;

public:
    Continuity_n(const Parameters &params);

    //!Sets up the matrix equation An*n = bn for continuity equation for electrons.
    //!\param V stores the voltage and is needed to calculate Bernoulli fnc.'s.
    //!\param Un stores the net generation rate, needed for the right hand side.
    //!\param n the electron density is needed to setup the boundary conditions.
    void setup_eqn(const Eigen::MatrixXd &V_matrix, const Eigen::MatrixXd &Un_matrix, const std::vector<double> &n);

    //setters for BC's:
    //for left and right BC's, will use input from the n matrix to determine
    void set_n_topBC();
    void set_n_bottomBC();
    void set_n_leftBC(const std::vector<double> &n);
    void set_n_rightBC(const std::vector<double> &n);

    //getters (const keyword ensures that fnc doesn't change anything)
    std::vector<double> get_main_diag() const {return main_diag;}
    //std::vector<double> get_upper_diag() const {return upper_diag;}
    //std::vector<double> get_lower_diag() const {return lower_diag;}
    //std::vector<double> get_rhs() const {return rhs;}
    //std::vector<std::vector<double> > get_n_mob() const {return n_mob;}
    Eigen::MatrixXd get_n_mob() const {return n_mob;}
    //std::vector<double> get_B_n1() const {return B_n1;}
    //std::vector<double> get_B_n2() const {return B_2;}
    Eigen::VectorXd get_rhs() const {return VecXd_rhs;}  //returns the Eigen object
    Eigen::SparseMatrix<double> get_sp_matrix() const {return sp_matrix;}

    std::vector<double> get_n_topBC() const {return n_topBC;}
    std::vector<double> get_n_bottomBC() const {return n_bottomBC;}
    std::vector<double> get_n_leftBC() const {return n_leftBC;}
    std::vector<double> get_n_rightBC() const {return n_rightBC;}

private:
    std::vector<double> far_lower_diag;
    std::vector<double> far_upper_diag;
    std::vector<double> main_diag;
    std::vector<double> upper_diag;
    std::vector<double> lower_diag;
    std::vector<double> rhs;
    Eigen::MatrixXd n_mob;  //!Matrix storing the position dependent electron mobility
    Eigen::VectorXd VecXd_rhs;  //rhs in Eigen object vector form, for sparse matrix solver
    Eigen::SparseMatrix<double> sp_matrix;

    //std::vector<Trp> triplet_list;
    //int trp_cnt;  //for counting the triplets

    //Boundary conditions
    std::vector<double> n_leftBC, n_rightBC, n_bottomBC, n_topBC;

    Eigen::MatrixXd Bn_posX;  //bernoulli (+dV_x)
    Eigen::MatrixXd Bn_negX;  //bernoulli (-dV_x)
    Eigen::MatrixXd Bn_posZ;  //bernoulli (+dV_z)
    Eigen::MatrixXd Bn_negZ;  //bernoulli (-dV_z)

    double Cn;
    int num_cell, num_elements;
    int N;

    //!Calculates the Bernoulli functions for dV in x direction and updates member arrays
    void Bernoulli_n_X(const Eigen::MatrixXd &V_matrix);

    //!Calculates the Bernoulli functions for dV in z direction and updates member arrays
    void Bernoulli_n_Z(const Eigen::MatrixXd &V_matrix);

    //matrix setup functions
    void set_far_lower_diag();
    void set_lower_diag();
    void set_main_diag();
    void set_upper_diag();
    void set_far_upper_diag();
    void set_rhs(const Eigen::MatrixXd &Un_matrix);
};

#endif // CONTINUITY_N_H
