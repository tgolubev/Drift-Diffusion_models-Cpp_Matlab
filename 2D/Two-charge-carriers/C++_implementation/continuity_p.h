#ifndef CONTINUITY_P_H
#define CONTINUITY_P_H

#include <vector>
#include <Eigen/Dense>

#include "parameters.h"  //needs this to know what parameters is
#include "constants.h"

class Continuity_p
{
public:   
    Continuity_p(const Parameters &params);

    //!Sets up the matrix equation Ap*p = bp for continuity equation for holes.
    //!\param V stores the voltage and is needed to calculate Bernoulli fnc.'s.
    //!\param Up stores the net generation rate, needed for the right hand side.
    void setup_eqn(const Eigen::MatrixXd &V_matrix, const Eigen::MatrixXd &Up_matrix);

    //setters for BC's:
    //for left and right BC's, will use input from the n matrix to determine
    std::vector<double> set_p_topBC();
    std::vector<double> set_p_bottomBC();
    std::vector<double> set_p_leftBC();
    std::vector<double> set_p_rightBC();

    //getters
    std::vector<double> get_main_diag() const {return main_diag;}
    std::vector<double> get_upper_diag() const {return upper_diag;}
    std::vector<double> get_lower_diag() const {return lower_diag;}
    std::vector<double> get_rhs() const {return rhs;}
    //std::vector<std::vector<double> > get_p_mob() const {return p_mob;}
    Eigen::MatrixXd get_p_mob() const {return p_mob;}

    //std::vector<double> get_B_p1() const {return B_p1;}
    //std::vector<double> get_B_p2() const {return B_p2;}

    std::vector<double> get_p_topBC() const {return p_topBC;}
    std::vector<double> get_p_bottomBC() const {return p_bottomBC;}
    std::vector<double> get_p_leftBC() const {return p_leftBC;}
    std::vector<double> get_p_rightBC() const {return p_rightBC;}

private:
    std::vector<double> far_lower_diag;
    std::vector<double> far_upper_diag;
    std::vector<double> main_diag;
    std::vector<double> upper_diag;
    std::vector<double> lower_diag;
    std::vector<double> rhs;
    Eigen::MatrixXd p_mob;  //!Matrix storing the position dependent holeelectron mobility

    //Boundary conditions
    std::vector<double> p_leftBC, p_rightBC, p_bottomBC, p_topBC;

    //Bernoulli functions
    Eigen::MatrixXd Bp_posX;  //bernoulli (+dV_x)
    Eigen::MatrixXd Bp_negX;  //bernoulli (-dV_x)
    Eigen::MatrixXd Bp_posZ;  //bernoulli (+dV_z)
    Eigen::MatrixXd Bp_negZ;  //bernoulli (-dV_z)

    double Cp;
    int num_cell, num_elements; //so don't have to keep typing params.
    int N;

    //!Calculates the Bernoulli functions for dV in x direction and updates member arrays
    void Bernoulli_p_X(const Eigen::MatrixXd &V_matrix);

    //!Calculates the Bernoulli functions for dV in z direction and updates member arrays
    void Bernoulli_p_Z(const Eigen::MatrixXd &V_matrix);

    //matrix setup functions
    void set_far_lower_diag();
    void set_lower_diag();
    void set_main_diag();
    void set_upper_diag();
    void set_far_upper_diag();
    void Continuity_p::set_rhs(const Eigen::MatrixXd &Up_matrix);
};

#endif // CONTINUITY_P_H
