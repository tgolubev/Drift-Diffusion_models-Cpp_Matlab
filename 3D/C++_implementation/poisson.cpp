#include "poisson.h"
#include <iostream>

Poisson::Poisson(const Parameters &params)
{
    CV = (params.N_dos*params.dx*params.dx*q)/(epsilon_0*Vt);
    N = params.num_cell -1;  //for convenience define this --> is the number of points in 1D inside the device
    num_elements = params.num_elements;
    num_cell = params.num_cell;
    V_matrix = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);    //useful for calculating currents at end of each Va
    V_matrix.setZero();
    netcharge = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
    netcharge.setZero();

    main_diag.resize(num_elements+1);
    upper_diag.resize(num_elements);
    lower_diag.resize(num_elements);
    far_lower_diag.resize(num_elements-N+1);
    far_upper_diag.resize(num_elements-N+1);
    rhs.resize(num_elements+1);  //+1 b/c I am filling from index 1

    V_bottomBC.resize(num_cell+1, num_cell+1);
    V_topBC.resize(num_cell+1, num_cell+1);
    V_leftBC_X.resize(num_cell+1, num_cell+1);
    V_rightBC_X.resize(num_cell+1, num_cell+1);
    V_leftBC_Y.resize(num_cell+1, num_cell+1);
    V_rightBC_Y.resize(num_cell+1, num_cell+1);

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

void Poisson::set_V_leftBC_X(const std::vector<double> &V)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            V_leftBC_X(j,k) = V[index + (j-1)*N + 1];
        }
        index = index+N*N;  //brings us to next vertical subblock set
    }
}

void Poisson::set_V_rightBC_X(const std::vector<double> &V)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            V_rightBC_X(j,k) = V[index + j*N];
        }
        index = index+N*N;  //brings us to next vertical subblock set
    }
}

void Poisson::set_V_leftBC_Y(const std::vector<double> &V)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int i = 1; i <= N; i++) {
            V_leftBC_Y(i, k) = V[index + i];
        }
        index = index + N*N;
    }
}

void Poisson::set_V_rightBC_Y(const std::vector<double> &V)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int i = 1; i <= N; i++) {
            V_rightBC_Y(i, k) = V[index + i + N*N - N];
        }
        index = index + N*N;
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

void Poisson::set_rhs(const std::vector<double> &n, const std::vector<double> &p)
{

    for (int i = 1; i <= num_elements; i++)
        rhs[i] = CV*(p[i] - n[i]);  //Note: this uses full device

    //add on BC's
    int index = 0;
    for (int k = 1; k <= N; k++) {
        if (k == 1) {
            for (int j = 1; j <= N; j++) {
                if (j == 1)  {//different for 1st subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)      //1st element has 2 BC's
                            rhs[index] += epsilon(i,j,k)*V_leftBC_X(1,1) + epsilon(i,j,k)*V_leftBC_Y(1,1) + epsilon(i,j,k)*V_bottomBC(i,j);
                        else if (i == N)
                            rhs[index] += epsilon(N+1, 1, 1)*V_rightBC_X(1,1) + epsilon(N, 0, 1)*V_leftBC_Y(N,1) + epsilon(N, 1, 0)*V_bottomBC(i,j);
                        else //middle elements in 1st subblock
                            rhs[index] += epsilon(i, 0, 1)*V_leftBC_Y(i, 1) + epsilon(i, 1, 0)*V_bottomBC(i,j);
                    }
                } else if (j == N) {      //different for last subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)  //1st element has 2 BC's
                            rhs[index] += epsilon(0, N, 1)*V_leftBC_X(N,1) + epsilon(1, N+1, 1)*V_rightBC_Y(1,1) + epsilon(1, N, 0)*V_bottomBC(i,j);
                        else if (i==N)
                            rhs[index] += epsilon(N+1, N, 1)*V_rightBC_X(N,1) + epsilon(N, N+1, 1)*V_rightBC_Y(N,1) + epsilon(N,N, 0)*V_bottomBC(i,j);
                        else //inner rows of Nth subblock
                            rhs[index] += epsilon(i,N+1, 1)*V_rightBC_Y(i,1) + epsilon(i, N, 0)*V_bottomBC(i,j);
                    }
                } else {     //interior subblocks of 1st (k = 1) subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)
                            rhs[index] += epsilon(0, j, 1)*V_leftBC_X(j,1) + epsilon(1, j, 0)*V_bottomBC(i,j);
                        else if (i == N)
                            rhs[index] += epsilon(N+1, j, 1)*V_rightBC_X(j,1) + epsilon(N, j, 0)*V_bottomBC(i,j);
                        else
                            rhs[index] += epsilon(i, j, 0)*V_bottomBC(i,j);
                    }
                }
            }
        } else if (k == N) {
            for (int j = 1; j <= N; j++) {
                if (j == 1)  {//different for 1st subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)      //1st element has 2 BC's
                            rhs[index] += epsilon(0,1,N)*V_leftBC_X(1,N) + epsilon(1,0,N)*V_leftBC_Y(1,N) + epsilon(1,1,N+1)*V_topBC(i,j);
                        else if (i == N)
                            rhs[index] += epsilon(N+1, 1, N)*V_rightBC_X(1,N) + epsilon(N, 0, N)*V_leftBC_Y(N,N) + epsilon(N, 1, N+1)*V_topBC(i,j);
                        else //middle elements in 1st subblock
                            rhs[index] += epsilon(i, 0, N)*V_leftBC_Y(i, N) + epsilon(i, 1, N+1)*V_topBC(i,j);
                    }
                } else if (j == N) {      //different for last subblock within k=N subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)  //1st element has 2 BC's
                            rhs[index] += epsilon(0,N, N)*V_leftBC_X(N,N) + epsilon(1, N+1, N)*V_rightBC_Y(1,N) + epsilon(1, N, N+1)*V_topBC(i,j);
                        else if (i==N)
                            rhs[index] += epsilon(N+1,N, N)*V_rightBC_X(N,N) + epsilon(N, N+1, N)*V_rightBC_Y(N,N) + epsilon(N,N, N+1)*V_topBC(i,j);
                        else //inner rows of Nth subblock
                            rhs[index] += epsilon(i,N+1, N)*V_rightBC_Y(i,N) + epsilon(i, N, N+1)*V_topBC(i,j);
                    }
                } else {    //interior subblocks of last (k = N) subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)
                            rhs[index] += epsilon(0, j, N)*V_leftBC_X(j,N) + epsilon(1, j, N+1)*V_topBC(i,j);
                        else if (i == N)
                            rhs[index] += epsilon(N+1, j, N)*V_rightBC_X(j,N) + epsilon(N, j, N+1)*V_topBC(i,j);
                        else
                            rhs[index] += epsilon(i, j, N+1)*V_topBC(i,j);
                    }
                }
            }
        } else {  //interior subblock groups (k=2:N-1)
            for (int j = 1; j <= N; j++) {
                if (j == 1)  {//different for 1st subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)      //1st element has 2 BC's
                            rhs[index] += epsilon(0, 1, k)*V_leftBC_X(1,k) + epsilon(1, 0, k)*V_leftBC_Y(1,k);
                        else if (i == N)
                            rhs[index] += epsilon(N+1, 1, k)*V_rightBC_X(1,k) + epsilon(N, 0, k)*V_leftBC_Y(N,k);
                        else //middle elements in 1st subblock
                            rhs[index] += epsilon(i, 0, k)*V_leftBC_Y(i, k);
                    }
                } else if (j == N) {      //different for last subblock within k=N subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)  //1st element has 2 BC's
                            rhs[index] += epsilon(0,N, k)*V_leftBC_X(N,k) + epsilon(1, N+1, k)*V_rightBC_Y(1,k);
                        else if (i==N)
                            rhs[index] += epsilon(N+1,N, k)*V_rightBC_X(N,k) + epsilon(N, N+1, k)*V_rightBC_Y(N,k);
                        else //inner rows of Nth subblock
                            rhs[index] += epsilon(i,N+1, k)*V_rightBC_Y(i,k);
                    }
                } else {    //interior subblocks of last (k = N) subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)
                            rhs[index] += epsilon(0, j, k)*V_leftBC_X(j,k);
                        else if (i == N)
                            rhs[index] += epsilon(N+1, j, k)*V_rightBC_X(j,k);
                    }
                }
            }
        }
    }

    //set up VectorXd Eigen vector object for sparse solver
    for (int i = 1; i<=num_elements; i++) {
        VecXd_rhs(i-1) = rhs[i];   //fill VectorXd  rhs of the equation
    }

}

//-----------------------------------------  //THIS NEEDS TO BE UPDATED-->
void Poisson::to_matrix(const std::vector<double> &V)
{
    for (int index = 1; index <= num_elements; index++) {
        int i = index % N;    //this gives the i value for matrix
        if (i == 0) i = N;

        int j = 1 + static_cast<int>(floor((index-1)/N));  // j value for matrix
        int k = 1 + static_cast<int>(floor((index-1)/(N*N)));

        V_matrix(i,j, k) = V[index];
    }

    for (int i = 0; i <= num_cell; i++) {
        for (int j = 0; j <= num_cell; j++) {
            V_matrix(i, j, 0) = V_bottomBC(i,j);
            V_matrix(i,j, num_cell) = V_topBC(i,j);
        }
    }

    for (int j = 1; j < num_cell; j++) {   //don't need to set j = 0 and j = num_cell elements, b/c already set when apply top and bottom BC's (are corners)
        for (int k = 1; k <= num_cell; k++) {
            V_matrix(0, j, k) = V_leftBC_X(j,k);
            V_matrix(num_cell, j, k) = V_rightBC_X(j,k);
        }
    }

    for (int i = 1; i < num_cell; i++) {   //don't need to set i = 0 and i = num_cell elements, b/c already set when apply top and bottom BC's (are corners)
        for (int k = 1; k <= num_cell; i++) {
            V_matrix(i, 0, k) = V_leftBC_Y(i,k);
            V_matrix(i, num_cell, k) = V_rightBC_Y(i,k);
        }
    }

}
