#ifndef THOMAS_TRIDIAG_SOLVE_H
#define THOMAS_TRIDIAG_SOLVE_H

#include<vector>

//!This function uses Thomas algorithm for tridiagonal matrix (special case of Gaussian elimination)
//! diagonal = array containing elements of main diagonal. indices: (a1.....an)
//! b = array containing elements of upper diagonal. indices (b1....b_n-1)
//!c = array containing elements of lower diagonal. indices (c1...c_n-1)
std::vector<double> Thomas_solve(const std::vector<double> &a, const std::vector<double> &b,const std::vector<double> &c, std::vector<double> &rhs);



#endif // THOMAS_TRIDIAG_SOLVE_H
