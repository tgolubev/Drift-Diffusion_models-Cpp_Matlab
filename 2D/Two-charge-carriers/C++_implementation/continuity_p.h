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
    void setup_eqn(const std::vector<double> &V, const std::vector<double> & Up);

    //getters
    std::vector<double> get_main_diag() const {return main_diag;}
    std::vector<double> get_upper_diag() const {return upper_diag;}
    std::vector<double> get_lower_diag() const {return lower_diag;}
    std::vector<double> get_rhs() const {return rhs;}
    //std::vector<std::vector<double> > get_p_mob() const {return p_mob;}
    Eigen::MatrixXd get_p_mob() const {return p_mob;}
    std::vector<double> get_B_p1() const {return B_p1;}
    std::vector<double> get_B_p2() const {return B_p2;}
    double get_p_bottomBC() const {return p_bottomBC;}
    double get_p_topBC() const {return p_topBC;}

private:
    std::vector<double> main_diag;
    std::vector<double> upper_diag;
    std::vector<double> lower_diag;
    std::vector<double> rhs;
    Eigen::MatrixXd p_mob;  //!Matrix storing the position dependent holeelectron mobility
    //std::vector<std::vector<double> > p_mob; //!Matrix storing the position dependent hole mobility
    std::vector<double> B_p1;  //bernoulli (+dV)
    std::vector<double> B_p2;  //bernoulli (-dV)
    double Cp;
    std::vector<double> p_leftBC, p_rightBC;
    double p_bottomBC, p_topBC;
    int num_cell; //so don't have to keep typing params.

    //!Calculates the Bernoulli functions and updates member arrays
    //! \param B_n1 = B(+dV) and \param B_n2 = (-dV)
    void BernoulliFnc_p(const std::vector<double> &V);

    void set_main_diag();
    void set_upper_diag();
    void set_lower_diag();
    void set_rhs(const std::vector<double> &Up);
};

#endif // CONTINUITY_P_H
