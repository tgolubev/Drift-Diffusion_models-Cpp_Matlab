#include "poisson.h"
#include <iostream>

Poisson::Poisson(const Parameters &params)
{
    CV = (params.N_dos*params.dx*params.dx*q)/(epsilon_0*Vt);
    Nx = params.Nx;
    Ny = params.Ny;
    Nz = params.Nz;
    num_elements = params.num_elements;
    V_matrix = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);    //useful for calculating currents at end of each Va
    V_matrix.setZero();
    netcharge = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
    netcharge.setZero();

    main_diag.resize(num_elements+1);
    upper_diag.resize(num_elements+1);
    lower_diag.resize(num_elements+1);
    far_lower_diag.resize(num_elements+1);
    far_upper_diag.resize(num_elements+1);
    rhs.resize(num_elements+1);  //+1 b/c I am filling from index 1

    V_bottomBC.resize(num_cell+1, num_cell+1);
    V_topBC.resize(num_cell+1, num_cell+1);

    //----------------------------------------------------------------------------------------------------------
    // //MUST FILL WITH THE VALUES OF epsilon!!  WILL NEED TO MODIFY THIS WHEN HAVE SPACE VARYING
    epsilon = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
    epsilon_avg_X = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
    epsilon_avg_Y = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
    epsilon_avg_Z = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
    epsilon.setConstant(params.eps_active/params.mobil);

    //Compute averaged mobilities
    for (int k = 1; k <= num_cell+1; k++) {
        for (int j = 1; j <= num_cell+1; j++) {
            for (int i = 1; i <= num_cell+1; i++) {
                epsilon_avg_X(i,j,k) = (epsilon(i,j,k) + epsilon(i,j+1,k) + epsilon(i,j,k+1) + epsilon(i,j+1,k+1))/4.;
                epsilon_avg_Y(i,j,k) = (epsilon(i,j,k) + epsilon(i+1,j,k) + epsilon(i,j,k+1) + epsilon(i+1,j,k+1))/4.;
                epsilon_avg_Z(i,j,k) = (epsilon(i,j,k) + epsilon(i+1,j,k) + epsilon(i,j+1,k) + epsilon(i+1,j+1,k))/4.;
            }
        }
    }
    //----------------------------------------------------------------------------------------------------------

    //allocate memory for the sparse matrix and rhs vector (Eig object)
    sp_matrix.resize(num_elements, num_elements);
    VecXd_rhs.resize(num_elements);   //only num_elements, b/c filling from index 0 (necessary for the sparse solver)

    //setup the triplet list for sparse matrix
    triplet_list.resize(7*num_elements);   //approximate the size that need

} //constructor)

//-------------------------------------------------------
//Set BC's functions


void Poisson::set_V_topBC(const Parameters &params, double Va)
{
    for (int j = 0; j <= num_cell; j++) {
        for (int i = 0; i <= num_cell; i++) {
            V_topBC(i,j) = (params.Vbi-Va)/(2*Vt) - params.phi_c/Vt;
        }
    }
}

void Poisson::set_V_bottomBC(const Parameters &params, double Va)
{
    for (int j = 0; j <= num_cell; j++) {
        for (int i = 0; i <= num_cell; i++) {
            V_bottomBC(i,j) = -((params.Vbi-Va)/(2*Vt) - params.phi_a/Vt);
        }
    }
}



void Poisson::setup_matrix()  //Note: this is on purpose different than the setup_eqn used for Continuity eqn's, b/c I need to setup matrix only once
{
    set_main_diag();
    set_lower_diag();
    set_upper_diag();
    set_far_lower_diag();
    set_far_upper_diag();

    typedef Eigen::Triplet<double> Trp;

    //generate triplets for Eigen sparse matrix
    //setup the triplet list for sparse matrix

     sp_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());    //sp_matrix is our sparse matrix
}


//---------------Setup AV diagonals (Poisson solve)---------------------------------------------------------------
void Poisson::set_far_lower_diag()
{
    int index = 1;
    for (int k = 2; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                triplet_list[trp_cnt] = {index-1+N*N, index-1, -epsilon_avg_Z(i,j,k)};  //note: don't need +1, b/c c++ values correspond directly to the inside pts
                //just  fill directly!! the triplet list. DON'T NEED THE DIAG VECTORS AT ALL!
                //RECALL, THAT the sparse matrices are indexed from 0 --> that's why have the -1's
                trp_cnt++;
                index = index +1;
            }
        }
    }
}

void Poisson::set_lower_diag()
{
    int index = 1;
    for (int k = 1; k <= N; k++) {
        for (int j = 2; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                triplet_list[trp_cnt] = {index-1+N, index-1, -epsilon_avg_Y(i,j,k)};
                trp_cnt++;
                index = index +1;
            }
        }
        index = index + N;  //add on the 0's subblock, so that filling it is skipped
    }
}


//main lower diag
void Poisson::set_main_lower_diag()
{
    int index = 1;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 2; i <= N; i++) {
                triplet_list[trp_cnt] = {index, index-1, -epsilon_avg_X(i,j,k)};
                trp_cnt++;
                index = index +1;
            }
            index = index + 1;  //skip the corner elements which are zero
        }
    }

}


void Poisson::set_main_diag()
{
    int index = 1;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                triplet_list[trp_cnt] = {index-1, index-1, epsilon_avg_X(i,j,k) + epsilon_avg_X(i+1,j,k) + epsilon_avg_Y(i,j,k) + epsilon_avg_Y(i,j+1,k) + epsilon_avg_Z(i,j,k) + epsilon_avg_Z(i,j,k+1)};
                trp_cnt++;
                index = index +1;
            }
        }
    }
}


void Poisson::set_main_upper_diag()
{

    int index = 1;  //note: unlike Matlab, can always start index at 1 here, b/c not using any spdiags fnc
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N-1; i++) {
                triplet_list[trp_cnt] = {index-1, index, -epsilon_avg_X(i+1,j,k)};
                trp_cnt++;
                index = index +1;
            }
            index = index +1;
        }
    }

}


void Poisson::set_upper_diag()
{

    int index = 1;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N-1; j++) {
            for (int i = 1; i <= N; i++) {
                triplet_list[trp_cnt] = {index-1, index-1+N, -epsilon_avg_Y(i,j+1,k)};
                trp_cnt++;
                index = index +1;
            }
        }
        index = index + N;
    }
}


void Poisson::set_far_upper_diag()
{

  int index = 1;
  for (int k = 1; k <= N-1; k++) {
      for (int j = 1; j <= N; j++) {
          for (int i = 1; i <= N; i++) {
               triplet_list[trp_cnt] = {index-1, index-1+N*N, -epsilon_avg_Z(i,j,k+1)};
               trp_cnt++;
               index = index +1;
          }
      }
  }
}

//---------------------------------------------------------------------------------------------------

void Poisson::set_rhs(const Eigen::Tensor<double, 3> &p)
{

    for (int i = 1; i <= num_elements; i++)
        rhs[i] = CV*(p[i]);  //Note: this uses full device

    //add on BC's


    //set up VectorXd Eigen vector object for sparse solver
    for (int i = 1; i<=num_elements; i++) {
        VecXd_rhs(i-1) = rhs[i];   //fill VectorXd  rhs of the equation
    }

}
