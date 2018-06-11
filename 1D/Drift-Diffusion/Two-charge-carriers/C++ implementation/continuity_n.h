#ifndef CONTINUITY_N_H
#define CONTINUITY_N_H

#include <vector>
#include "parameters.h"  //needs this to know what parameters are
#include "constants.h"

class Continuity_n
{
public:
    Continuity_n(Parameters &params);
    void setup_eqn(std::vector<double> &B_n1, std::vector<double> &B_n2, std::vector<double> &Un);

    //getters
    std::vector<double> get_main_diag() const {return main_diag;}  //const keyword ensures that fnc doesn't change anything
    std::vector<double> get_upper_diag() const {return upper_diag;}
    std::vector<double> get_lower_diag() const {return lower_diag;}
    std::vector<double> get_rhs() const {return rhs;}
    std::vector<double> get_n_mob() const {return n_mob;}
    double get_n_leftBC() const {return n_leftBC;}
    double get_n_rightBC() const {return n_rightBC;}

private:
    std::vector<double> main_diag;
    std::vector<double> upper_diag;
    std::vector<double> lower_diag;
    std::vector<double> rhs;
    std::vector<double> n_mob;
    double Cn;
    double n_leftBC;        //this is anode
    double n_rightBC;

    void set_main_diag(const  std::vector<double> &B_n1,const  std::vector<double> &B_n2);
    void set_upper_diag(const  std::vector<double> &B_n1);
    void set_lower_diag(const  std::vector<double> &B_n2);
    void set_rhs(const std::vector<double> &B_n1, const std::vector<double> &B_n2, std::vector<double> &Un);
};

#endif // CONTINUITY_N_H
