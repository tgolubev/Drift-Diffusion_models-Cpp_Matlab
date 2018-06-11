#ifndef CONTINUITY_P_H
#define CONTINUITY_P_H

#include <vector>
#include "parameters.h"  //needs this to know what parameters is
#include "constants.h"

class Continuity_p
{
public:   
    Continuity_p(Parameters &params);
    void setup_eqn(std::vector<double> &B_p1, std::vector<double> & B_p2, std::vector<double> & Up);

    //getters
    std::vector<double> get_main_diag() const {return main_diag;}
    std::vector<double> get_upper_diag() const {return upper_diag;}
    std::vector<double> get_lower_diag() const {return lower_diag;}
    std::vector<double> get_rhs() const {return rhs;}
    std::vector<double> get_p_mob() const {return p_mob;}
    double get_p_leftBC() const {return p_leftBC;}
    double get_p_rightBC() const {return p_rightBC;}

private:
    std::vector<double> main_diag;
    std::vector<double> upper_diag;
    std::vector<double> lower_diag;
    std::vector<double> rhs;
    std::vector<double> p_mob;
    double Cp;
    double p_leftBC;
    double p_rightBC;

    void set_main_diag(const  std::vector<double> &B_p1,const  std::vector<double> &B_p2);
    void set_upper_diag(const  std::vector<double> &B_p1);
    void set_lower_diag(const  std::vector<double> &B_p2);
    void set_rhs(const  std::vector<double> &B_p1, const std::vector<double> &B_p2, std::vector<double> &Up);
};

#endif // CONTINUITY_P_H
