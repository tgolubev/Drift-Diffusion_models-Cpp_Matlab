#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include "parameters.h"
#include <fstream>
#include <iomanip>
#include <iostream>

#include <Eigen/Dense>

class Utilities
{
public:
    Utilities();

    //!Applies linear mixing (result = w*new_value + (1-w)*old_value) where w is the mixing factor.
    //! The mixing factor is in the \param params object.
    std::vector<double> linear_mix(const Parameters &params,const std::vector<double> &new_values,const std::vector<double> &old_values);

    //!This writes to output files the details of voltage \param V, carrier densities \param p and \param n, current \param J_total, net electron generation rate \param Un.
    //! The files are named according to the applied voltage \param Va of this data.
    void write_details(const Parameters &params, double Va, const Eigen::MatrixXd &V_matrix, const Eigen::MatrixXd &p_matrix, const Eigen::MatrixXd &n_matrix, const Eigen::MatrixXd &J_total_Z, const Eigen::MatrixXd  &Un_matrix);

    //!Writes to a file the JV data. The file contains 3 columns, the applied voltage \param Va, the current \param J_total,
    //! and the number of iterations \param iter required to converge at that voltage.
    void write_JV(const Parameters &params, std::ofstream &JV, double iter, double Va, const Eigen::MatrixXd &J_total_Z);

};

#endif // UTILITIES_H
