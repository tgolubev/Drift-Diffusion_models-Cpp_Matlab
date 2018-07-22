#ifndef CONTINUITY_P_H
#define CONTINUITY_P_H

#include <vector>
#include <Eigen/Dense>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include<Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include <unsupported/Eigen/CXX11/Tensor>

#include "parameters.h"  //needs this to know what parameters is
#include "constants.h"

class Continuity_p
{

typedef Eigen::Triplet<double> Trp;  //allows to use Trp to refer to the Eigen::Triplet<double> type

public:
    Continuity_p(const Parameters &params);

    //!Sets up the matrix equation Ap*p = bp for continuity equation for holes.
    //!\param V stores the voltage and is needed to calculate Bernoulli fnc.'s.
    //!\param Up stores the net generation rate, needed for the right hand side.
    //!\param p the hole density is needed to setup the boundary conditions.
    void setup_eqn(const Eigen::Tensor<double, 3> &V_matrix, const std::vector<double> &Up, const Eigen::Tensor<double, 3> &p);

    void calculate_currents();

    //setters for BC's:
    //for left and right BC's, will use input from the n matrix to determine
    void set_p_topBC();
    void set_p_bottomBC();
    void set_p_leftBC_X(const std::vector<double> &p);
    void set_p_rightBC_X(const std::vector<double> &p);
    void set_p_leftBC_Y(const std::vector<double> &p);
    void set_p_rightBC_Y(const std::vector<double> &p);


    //getters
    Eigen::VectorXd get_rhs() const {return VecXd_rhs;}  //returns the Eigen object
    Eigen::SparseMatrix<double> get_sp_matrix() const {return sp_matrix;}

    double get_p_bottomBC(int i, int j) const {return p_bottomBC(i,j);}  //bottom and top are needed to set initial conditions
    double get_p_topBC(int i, int j) const {return p_topBC(i,j);}

    Eigen::Tensor<double, 3> get_p_matrix() const {return p_matrix;}
    Eigen::Tensor<double, 3> get_Jp_X() const {return Jp_X;}
    Eigen::Tensor<double, 3> get_Jp_Y() const {return Jp_Y;}
    Eigen::Tensor<double, 3> get_Jp_Z() const {return Jp_Z;}

    //The below getters can be useful for testing and debugging
    //std::vector<double> get_main_diag() const {return main_diag;}
    //std::vector<double> get_upper_diag() const {return upper_diag;}
    //std::vector<double> get_lower_diag() const {return lower_diag;}
    //std::vector<double> get_far_upper_diag() const {return far_upper_diag;}
    //std::vector<double> get_far_lower_diag() const {return far_lower_diag;}
    //Eigen::MatrixXd  get_Bp_posX() const {return Bp_posX;}
    //Eigen::MatrixXd  get_Bp_negX() const {return Bp_negX;}
    //Eigen::MatrixXd  get_Bp_posZ() const {return Bp_posZ;}
    //Eigen::MatrixXd  get_Bp_negZ() const {return Bp_negZ;}
    //Eigen::MatrixXd get_p_mob() const {return p_mob;}


private:
    //vectors for the 11 diagonals
     //I think I can fill the triplet list directly, and don't need these diag vectors
//    std::vector<double> lowest_diag;
//    std::vector<double> lower_diag_Xs;
//    std::vector<double> lower_diag_Y_PBCs;
//    std::vector<double> lower_diag_Ys;
//    std::vector<double> main_lower_diag;
//    std::vector<double> main_diag;
//    std::vector<double> main_upper_diag;
//    std::vector<double> upper_diag_Ys;
//    std::vector<double> upper_diag_Y_PBCs;
//    std::vector<double> upper_diag_Xs;
//    std::vector<double> highest_diag;

    std::vector<double> rhs;

    Eigen::Tensor<double, 3> p_mob;  //!Matrix storing the position dependent holeelectron mobility
    Eigen::Tensor<double, 3> p_mob_avg_X, p_mob_avg_Y, p_mob_avg_Z;
    Eigen::VectorXd VecXd_rhs;  //rhs in Eigen object vector form, for sparse matrix solver
    Eigen::SparseMatrix<double> sp_matrix;
    Eigen::Tensor<double, 3> p_matrix;
    Eigen::Tensor<double, 3> Jp_Z;
    Eigen::Tensor<double, 3> Jp_X;
    Eigen::Tensor<double, 3> Jp_Y;

    std::vector<Trp> triplet_list;
    int trp_cnt;  //for counting the triplets
    double J_coeff;  //coefficient for curents eqn

    Eigen::Tensor<double, 3> values;  //temporary variable to store values for diagonal filling
    Eigen::Tensor<double, 3> values2, values3, values4;  //these additional temp variables needed for main diag filling

    //Boundary conditions
    Eigen::MatrixXd p_leftBC_X, p_rightBC_X, p_leftBC_Y, p_rightBC_Y, p_bottomBC, p_topBC;

    //Bernoulli functions
    Eigen::Tensor<double, 3> Bp_posX;  //bernoulli (+dV_x)
    Eigen::Tensor<double, 3> Bp_negX;  //bernoulli (-dV_x)
    Eigen::Tensor<double, 3> Bp_posY;  //bernoulli (+dV_y)
    Eigen::Tensor<double, 3> Bp_negY;  //bernoulli (-dV_y)
    Eigen::Tensor<double, 3> Bp_posZ;  //bernoulli (+dV_z)
    Eigen::Tensor<double, 3> Bp_negZ;  //bernoulli (-dV_z)

    double Cp;
    int num_elements; //so don't have to keep typing params.
    int Nx, Ny, Nz;
    int num_cell_x, num_cell_y, num_cell_z;

    //!Calculates the Bernoulli functions for dVs and updates member arrays
    void Bernoulli_p(const Eigen::Tensor<double, 3> &V_matrix);

    //functions for setting up the 11 diagonals
    void set_lowest_diag();
    void set_lower_diag_Xs();   //lower diag corresponding to X direction finite differences
    void set_lower_diag_Y_PBCs(); //lower diag corresponding to Y periodic boundary conditions
    void set_lower_diag_Ys();
    void set_main_lower_diag();
    void set_main_diag();
    void set_main_upper_diag();  //corresponds to Z direction finite differences
    void set_upper_diag_Ys();
    void set_upper_diag_Y_PBCs();
    void set_upper_diag_Xs();
    void set_highest_diag();

    void Continuity_p::set_rhs(const std::vector<double> &Up);
};

#endif // CONTINUITY_P_H
