#ifndef SIMULATION_H
#define SIMULATION_H


class Simulation
{
public:
    Simulation();
    void reduce_w(){w = w/w_reduce_factor;}
    void relax_tolerance(){tolerance = tolerance*tol_relax_factor;}
    void use_tolerance_eq() {tolerance = tolerance_eq;}
    void use_tolerance_i() {tolerance = tolerance_i;}
    void use_w_i() {w = w_i;}
    void use_w_eq() {w = w_eq;}

    //getters
    double get_tolerance() const {return tolerance;}
    double get_w() const {return w;}   //const specifies that fnc shouldn't change w
    int get_num_cell() const {return num_cell;}
    int get_num_V() const {return num_V;}
    double get_Va_min() const {return Va_min;}
    double get_Va_max() const {return Va_max;}
    double get_increment() const {return increment;}
    //double dx() const {return dx;}

private:
    double w;
    double w_reduce_factor;
    double tolerance;
    double tol_relax_factor;
    double Vmin, Vmax;
    double tolerance_eq, tolerance_i, w_i, w_eq;
    double L;
    int num_cell;
    //const double dx;

    double Va_min, Va_max, increment;
    int num_V;

};

#endif // SIMULATION_H
