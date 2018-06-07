#ifndef CONTINUITY_N_H
#define CONTINUITY_N_H

#include <vector>
#include "parameters.h"  //needs this to know what parameters are
#include "constants.h"

class Continuity_n
{
public:
    Continuity_n(Parameters &params): main_diag(params.get_num_cell()), upper_diag(params.get_num_cell()-1), lower_diag(params.get_num_cell()-1),  rhs(params.get_num_cell())
    {
        Cn = params.dx*params.dx/(Vt*params.N*params.mobil);
        n_leftBC = (params.N_LUMO*exp(-(params.E_gap - params.phi_a)/Vt))/params.N;       //this is anode
        n_rightBC = (params.N_LUMO*exp(-params.phi_c/Vt))/params.N;
    }

    void setup_eqn(std::vector<double> &n_mob, std::vector<double> &B_n1, std::vector<double> &B_n2, std::vector<double> &Un);

    //getters
    std::vector<double> get_main_diag() const {return main_diag;}  //const keyword ensures that fnc doesn't change anything
    std::vector<double> get_upper_diag() const {return upper_diag;}
    std::vector<double> get_lower_diag() const {return lower_diag;}
    std::vector<double> get_rhs() const {return rhs;}
    double get_n_leftBC() const {return n_leftBC;}
    double get_n_rightBC() const {return n_rightBC;}

private:
    std::vector<double> main_diag;
    std::vector<double> upper_diag;
    std::vector<double> lower_diag;
    std::vector<double> rhs;
    double Cn;
    double n_leftBC;        //this is anode
    double n_rightBC;

    void set_main_diag(const std::vector<double> &n_mob,const  std::vector<double> &B_n1,const  std::vector<double> &B_n2);
    void set_upper_diag(const std::vector<double> &n_mob,const  std::vector<double> &B_n1);
    void set_lower_diag(const std::vector<double> &n_mob,const  std::vector<double> &B_n2);
    void set_rhs(const std::vector<double> &n_mob, const std::vector<double> &B_n1, const std::vector<double> &B_n2, std::vector<double> &Un);
};

#endif // CONTINUITY_N_H
