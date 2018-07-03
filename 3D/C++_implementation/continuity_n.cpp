#include "continuity_n.h"

Continuity_n::Continuity_n(const Parameters &params)
{
    num_elements = params.num_elements; //note: num_elements is same thing as num_rows in main.cpp
    N = params.num_cell - 1;
    num_cell = params.num_cell;
    n_matrix = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);    //seems this is the only way to resize a tensor...

    main_diag.resize(num_elements+1);
    main_upper_diag.resize(num_elements+1);
    main_lower_diag.resize(num_elements+1);
    upper_diag.resize(num_elements+1);
    lower_diag.resize(num_elements+1);
    far_lower_diag.resize(num_elements+1);
    far_upper_diag.resize(num_elements+1);
    rhs.resize(num_elements+1);  //+1 b/c I am filling from index 1

   n_bottomBC.resize(num_cell+1, num_cell+1);
   n_topBC.resize(num_cell+1, num_cell+1);
   n_leftBC_X.resize(num_cell+1, num_cell+1);
   n_rightBC_X.resize(num_cell+1, num_cell+1);
   n_leftBC_Y.resize(num_cell+1, num_cell+1);
   n_rightBC_Y.resize(num_cell+1, num_cell+1);

   Bn_posX = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
   Bn_negX = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
   Bn_posY = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
   Bn_negY = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
   Bn_posZ = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
   Bn_negZ = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);

   values = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
   values2 = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
   values3 = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
   values4 = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);

   Jn_Z = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
   Jn_X = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
   Jn_Y = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);

   J_coeff = (q*Vt*params.N_dos*params.mobil)/params.dx;

   //------------------------------------------------------------------------------------------
   // //MUST FILL WITH THE VALUES OF n_mob!!  WILL NEED TO MODIFY THIS WHEN HAVE SPACE VARYING
   n_mob = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
   n_mob_avg_X = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
   n_mob_avg_Y = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
   n_mob_avg_Z = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
   n_mob.setConstant(params.n_mob_active/params.mobil);

   //Compute averaged mobilities
   for (int k = 1; k <= num_cell+1; k++) {
       for (int j = 1; j <= num_cell+1; j++) {
           for (int i = 1; i <= num_cell+1; i++) {
               n_mob_avg_X(i,j,k) = (n_mob(i,j,k) + n_mob(i,j+1,k) + n_mob(i,j,k+1) + n_mob(i,j+1,k+1))/4.;
               n_mob_avg_Y(i,j,k) = (n_mob(i,j,k) + n_mob(i+1,j,k) + n_mob(i,j,k+1) + n_mob(i+1,j,k+1))/4.;
               n_mob_avg_Z(i,j,k) = (n_mob(i,j,k) + n_mob(i+1,j,k) + n_mob(i,j+1,k) + n_mob(i+1,j+1,k))/4.;
           }
       }
   }
   //------------------------------------------------------------------------------------------
   Cn = (params.dx*params.dx)/(Vt*params.N_dos*params.mobil);

   //these BC's for now stay constant throughout simulation, so fill them once, upon Continuity_n object construction
   for (int j =  0; j <= num_cell; j++) {
       for (int i = 1; i <= num_cell+1; i++) {
           n_bottomBC(i,j) = params.N_LUMO*exp(-(params.E_gap-params.phi_a)/Vt)/params.N_dos;
           n_topBC(i,j) = params.N_LUMO*exp(-params.phi_c/Vt)/params.N_dos;
       }
   }

   //allocate memory for the sparse matrix and rhs vector (Eig object)
   sp_matrix.resize(num_elements, num_elements);
   VecXd_rhs.resize(num_elements);   //only num_elements, b/c filling from index 0 (necessary for the sparse solver)

   //setup the triplet list for sparse matrix
    triplet_list.resize(7*num_elements);   //approximate the size that need         // list of non-zeros coefficients in triplet form(row index, column index, value)
}

//----------------------------------------------------------
//Set BC's

void Continuity_n::set_n_leftBC_X(const std::vector<double> &n)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            n_leftBC_X(j,k) = n[index + (j-1)*N + 1];
        }
        index = index+N*N;  //brings us to next vertical subblock set
    }
}

void Continuity_n::set_n_rightBC_X(const std::vector<double> &n)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
         n_rightBC_X(j,k) = n[index + j*N];
        }
        index = index+N*N;  //brings us to next vertical subblock set
    }
}

void Continuity_n::set_n_leftBC_Y(const std::vector<double> &n)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int i = 1; i <= N; i++) {
            n_leftBC_Y(i, k) = n[index + i];
        }
        index = index + N*N;
    }
}

void Continuity_n::set_n_rightBC_Y(const std::vector<double> &n)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int i = 1; i <= N; i++) {
            n_rightBC_Y(i, k) = n[index + i + N*N - N];
        }
        index = index + N*N;
    }
}


//Calculates Bernoulli fnc values, then sets the diagonals and rhs
//use the V_matrix for setup, to be able to write equations in terms of (x,z) coordingates
void Continuity_n::setup_eqn(const Eigen::Tensor<double, 3> &V_matrix, const std::vector<double> &Un, const std::vector<double> &n)
{
    trp_cnt = 0;  //reset triplet count
    Bernoulli_n_X(V_matrix);
    Bernoulli_n_Y(V_matrix);
    Bernoulli_n_Z(V_matrix);

    set_far_lower_diag();
    set_lower_diag();
    set_main_lower_diag();
    set_main_diag();
    set_main_upper_diag();
    set_upper_diag();
    set_far_upper_diag();

    set_n_rightBC_X(n);
    set_n_leftBC_X(n);
    set_n_rightBC_Y(n);
    set_n_leftBC_Y(n);

    set_rhs(Un);

    sp_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());   //sp_matrix is our sparse matrix

}

//-------------------------------Setup An diagonals (Continuity/drift-diffusion solve)-----------------------------

void Continuity_n::set_far_lower_diag()
{
    values = -n_mob_avg_Z*Bn_negZ;  //element wise multiplication

    int index = 1;
    for (int k = 2; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                triplet_list[trp_cnt] = {index-1+N*N, index-1, values(i,j,k)};  //note: don't need +1, b/c c++ values correspond directly to the inside pts
                //just  fill directly!! the triplet list. DON'T NEED THE DIAG VECTORS AT ALL!
                //RECALL, THAT the sparse matrices are indexed from 0 --> that's why have the -1's
                trp_cnt++;
                index = index +1;
            }
        }
    }
}

void Continuity_n::set_lower_diag()
{
    values = -n_mob_avg_Y*Bn_negY;

    int index = 1;
    for (int k = 1; k <= N; k++) {
        for (int j = 2; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                triplet_list[trp_cnt] = {index-1+N, index-1, values(i,j,k)};
                trp_cnt++;
                index = index +1;
            }
        }
        index = index + N;  //add on the 0's subblock, so that filling it is skipped
    }
}


//main lower diag
void Continuity_n::set_main_lower_diag()
{
    values = -n_mob_avg_X*Bn_negX;

    int index = 1;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 2; i <= N; i++) {
                triplet_list[trp_cnt] = {index, index-1, values(i,j,k)};
                trp_cnt++;
                index = index +1;
            }
            index = index + 1;  //skip the corner elements which are zero
        }
    }

}


void Continuity_n::set_main_diag()
{
    values  = n_mob_avg_Z*Bn_posZ + n_mob_avg_Y*Bn_posY + n_mob_avg_X*Bn_posX;
    values2 = n_mob_avg_X*Bn_negX;
    values3 = n_mob_avg_Y*Bn_negY;
    values4 = n_mob_avg_Z*Bn_negZ;

    int index = 1;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                triplet_list[trp_cnt] = {index-1, index-1, values(i,j,k) + values2(i+1,j,k) + values3(i,j+1,k) + values4(i,j,k+1)};
                trp_cnt++;
                index = index +1;
            }
        }
    }
}


void Continuity_n::set_main_upper_diag()
{
    values = -n_mob_avg_X*Bn_posX;

    int index = 1;  //note: unlike Matlab, can always start index at 1 here, b/c not using any spdiags fnc
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N-1; i++) {
                triplet_list[trp_cnt] = {index-1, index, values(i+1,j,k)};
                trp_cnt++;
                index = index +1;
            }
            index = index +1;
        }
    }

}


void Continuity_n::set_upper_diag()
{
    values = -n_mob_avg_Y*Bn_posY;

    int index = 1;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N-1; j++) {
            for (int i = 1; i <= N; i++) {
                triplet_list[trp_cnt] = {index-1, index-1+N, values(i,j+1,k)};
                trp_cnt++;
                index = index +1;
            }
        }
        index = index + N;
    }
}


void Continuity_n::set_far_upper_diag()
{
  values = -n_mob_avg_Z*Bn_posZ;

  int index = 1;
  for (int k = 1; k <= N-1; k++) {
      for (int j = 1; j <= N; j++) {
          for (int i = 1; i <= N; i++) {
               triplet_list[trp_cnt] = {index-1, index-1+N*N, values(i,j,k+1)};
               trp_cnt++;
               index = index +1;
          }
      }
  }


}

//---------------------------------------------------------------------------

void Continuity_n::set_rhs(const std::vector<double> &Un)
{
    //calculate main part here
    for (int i = 1; i <= num_elements; i++)
        rhs[i] = Cn*Un[i];

    //add on BC's
    int index = 0;
    for (int k = 1; k <= N; k++) {
        if (k == 1) {
            for (int j = 1; j <= N; j++) {
                if (j == 1)  {//different for 1st subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)      //1st element has 2 BC's
                            rhs[index] += n_mob(i,j,k)*n_leftBC_X(1,1)*Bn_negX(i,j,k) + n_mob(i,j,k)*n_leftBC_Y(1,1)*Bn_negY(i,j,k) + n_mob(i,j,k)*n_bottomBC(i,j)*Bn_negZ(i,j,k);
                        else if (i == N)
                            rhs[index] += n_mob(N+1, 1, 1)*n_rightBC_X(1,1)*Bn_posX(i+1,j,k) + n_mob(N, 0, 1)*n_leftBC_Y(N,1)*Bn_negY(i,j,k) + n_mob(N, 1, 0)*n_bottomBC(i,j)*Bn_negZ(i,j,k);
                        else //middle elements in 1st subblock
                            rhs[index] += n_mob(i, 0, 1)*n_leftBC_Y(i, 1)*Bn_negY(i,j,k) + n_mob(i, 1, 0)*n_bottomBC(i,j)*Bn_negZ(i,j,k);
                    }
                } else if (j == N) {      //different for last subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)  //1st element has 2 BC's
                            rhs[index] += n_mob(0, N, 1)*n_leftBC_X(N,1)*Bn_negX(i,j,k) + n_mob(1, N+1, 1)*n_rightBC_Y(1,1)*Bn_posY(i,j+1,k) + n_mob(1, N, 0)*n_bottomBC(i,j)*Bn_negZ(i,j,k);
                        else if (i==N)
                            rhs[index] += n_mob(N+1, N, 1)*n_rightBC_X(N,1)*Bn_posX(i+1,j,k) + n_mob(N, N+1, 1)*n_rightBC_Y(N,1)*Bn_posY(i,j+1,k) + n_mob(N,N, 0)*n_bottomBC(i,j)*Bn_negZ(i,j,k);
                        else //inner rows of Nth subblock
                            rhs[index] += n_mob(i,N+1, 1)*n_rightBC_Y(i,1)*Bn_posY(i,j+1,k) + n_mob(i, N, 0)*n_bottomBC(i,j)*Bn_negZ(i,j,k);
                    }
                } else {     //interior subblocks of 1st (k = 1) subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)
                            rhs[index] += n_mob(0, j, 1)*n_leftBC_X(j,1)*Bn_negX(i,j,k) + n_mob(1, j, 0)*n_bottomBC(i,j)*Bn_negZ(i,j,k);
                        else if (i == N)
                            rhs[index] += n_mob(N+1, j, 1)*n_rightBC_X(j,1)*Bn_posX(i+1,j,k) + n_mob(N, j, 0)*n_bottomBC(i,j)*Bn_negZ(i,j,k);
                        else
                            rhs[index] += n_mob(i, j, 0)*n_bottomBC(i,j)*Bn_negZ(i,j,k);
                    }
                }
            }
        } else if (k == N) {
            for (int j = 1; j <= N; j++) {
                if (j == 1)  {//different for 1st subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)      //1st element has 2 BC's
                            rhs[index] += n_mob(0,1,N)*n_leftBC_X(1,N)*Bn_negX(i,j,k) + n_mob(1,0,N)*n_leftBC_Y(1,N)*Bn_negY(i,j,k) + n_mob(1,1,N+1)*n_topBC(i,j)*Bn_posZ(i,j,k+1);
                        else if (i == N)
                            rhs[index] += n_mob(N+1, 1, N)*n_rightBC_X(1,N)*Bn_posX(i+1,j,k) + n_mob(N, 0, N)*n_leftBC_Y(N,N)*Bn_negY(i,j,k) + n_mob(N, 1, N+1)*n_topBC(i,j)*Bn_posZ(i,j,k+1);
                        else //middle elements in 1st subblock
                            rhs[index] += n_mob(i, 0, N)*n_leftBC_Y(i, N)*Bn_negY(i,j,k) + n_mob(i, 1, N+1)*n_topBC(i,j)*Bn_posZ(i,j,k+1);
                    }
                } else if (j == N) {      //different for last subblock within k=N subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)  //1st element has 2 BC's
                            rhs[index] += n_mob(0,N, N)*n_leftBC_X(N,N)*Bn_negX(i,j,k) + n_mob(1, N+1, N)*n_rightBC_Y(1,N)*Bn_posY(i,j+1,k) + n_mob(1, N, N+1)*n_topBC(i,j)*Bn_posZ(i,j,k+1);
                        else if (i==N)
                            rhs[index] += n_mob(N+1,N, N)*n_rightBC_X(N,N)*Bn_posX(i+1,j,k) + n_mob(N, N+1, N)*n_rightBC_Y(N,N)*Bn_posY(i,j+1,k) + n_mob(N,N, N+1)*n_topBC(i,j)*Bn_posZ(i,j,k+1);
                        else //inner rows of Nth subblock
                            rhs[index] += n_mob(i,N+1, N)*n_rightBC_Y(i,N)*Bn_posY(i,j+1,k) + n_mob(i, N, N+1)*n_topBC(i,j)*Bn_posZ(i,j,k+1);
                    }
                } else {    //interior subblocks of last (k = N) subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)
                            rhs[index] += n_mob(0, j, N)*n_leftBC_X(j,N)*Bn_negX(i,j,k) + n_mob(1, j, N+1)*n_topBC(i,j)*Bn_posZ(i,j,k+1);
                        else if (i == N)
                            rhs[index] += n_mob(N+1, j, N)*n_rightBC_X(j,N)*Bn_posX(i+1,j,k) + n_mob(N, j, N+1)*n_topBC(i,j)*Bn_posZ(i,j,k+1);
                        else
                            rhs[index] += n_mob(i, j, N+1)*n_topBC(i,j)*Bn_posZ(i,j,k+1);
                    }
                }
            }
        } else {  //interior subblock groups (k=2:N-1)
            for (int j = 1; j <= N; j++) {
                if (j == 1)  {//different for 1st subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)      //1st element has 2 BC's
                            rhs[index] += n_mob(0, 1, k)*n_leftBC_X(1,k)*Bn_negX(i,j,k)  + n_mob(1, 0, k)*n_leftBC_Y(1,k)*Bn_negY(i,j,k);
                        else if (i == N)
                            rhs[index] += n_mob(N+1, 1, k)*n_rightBC_X(1,k)*Bn_posX(i+1,j,k) + n_mob(N, 0, k)*n_leftBC_Y(N,k)*Bn_negY(i,j,k);
                        else //middle elements in 1st subblock
                            rhs[index] += n_mob(i, 0, k)*n_leftBC_Y(i, k)*Bn_negY(i,j,k);
                    }
                } else if (j == N) {      //different for last subblock within k=N subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)  //1st element has 2 BC's
                            rhs[index] += n_mob(0,N, k)*n_leftBC_X(N,k)*Bn_negX(i,j,k) + n_mob(1, N+1, k)*n_rightBC_Y(1,k)*Bn_posY(i,j+1,k);
                        else if (i==N)
                            rhs[index] += n_mob(N+1,N, k)*n_rightBC_X(N,k)*Bn_posX(i+1,j,k) + n_mob(N, N+1, k)*n_rightBC_Y(N,k)*Bn_posY(i,j+1,k);
                        else //inner rows of Nth subblock
                            rhs[index] += n_mob(i,N+1, k)*n_rightBC_Y(i,k)*Bn_posY(i,j+1,k);
                    }
                } else {    //interior subblocks of last (k = N) subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)
                            rhs[index] += n_mob(0, j, k)*n_leftBC_X(j,k)*Bn_negX(i,j,k) ;
                        else if (i == N)
                            rhs[index] += n_mob(N+1, j, k)*n_rightBC_X(j,k)*Bn_posX(i+1,j,k);
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

//------------------------
//Note: are using the V matrix for Bernoulli calculations.
//Makes it clearer to write indices in terms of (x,z) real coordinate values.

void Continuity_n::Bernoulli_n_X(const Eigen::Tensor<double, 3> &V_matrix)
{
    Eigen::Tensor<double, 3> dV(num_cell+1, num_cell+1, num_cell+1);
    dV.setZero();

    for (int i = 1; i < num_cell+1; i++)
       for (int j = 1; j < num_cell+1; j++)
           for (int k = 1; k < num_cell+1; k++)
                dV(i,j,k) =  V_matrix(i,j,k)-V_matrix(i-1,j,k);

    for (int i = 1; i < num_cell+1; i++) {
        for (int j = 1; j < num_cell+1; j++) {
            for (int k = 1; k < num_cell+1; k++) {
                if (abs(dV(i,j,k)) < 1e-13) {        //to prevent blowup due  to 0 denominator
                    Bn_posX(i,j,k) = 1;//1 - dV(i,j)/2. + (dV(i,j)*dV(i,j))/12. - pow(dV(i,j), 4)/720.;
                    Bn_negX(i,j,k) =  1;//Bn_posX(i,j)*exp(dV(i,j));
                } else {
                    Bn_posX(i,j,k) = dV(i,j,k)/(exp(dV(i,j,k)) - 1.0);
                    Bn_negX(i,j,k) = Bn_posX(i,j,k)*exp(dV(i,j,k));
                }
            }
        }
    }
}

void Continuity_n::Bernoulli_n_Y(const Eigen::Tensor<double, 3> &V_matrix)
{
    Eigen::Tensor<double, 3> dV(num_cell+1, num_cell+1, num_cell+1);
    dV.setZero();

    for (int i = 1; i < num_cell+1; i++)
       for (int j = 1; j < num_cell+1; j++)
           for (int k = 1; k < num_cell+1; k++)
                dV(i,j,k) =  V_matrix(i,j,k)-V_matrix(i,j-1,k);

    for (int i = 1; i < num_cell+1; i++) {
        for (int j = 1; j < num_cell+1; j++) {
            for (int k = 1; k < num_cell+1; k++) {
                if (abs(dV(i,j,k)) < 1e-13) {        //to prevent blowup due  to 0 denominator
                    Bn_posY(i,j,k) = 1;//1 - dV(i,j)/2. + (dV(i,j)*dV(i,j))/12. - pow(dV(i,j), 4)/720.;
                    Bn_negY(i,j,k) =  1;//Bn_posZ(i,j)*exp(dV(i,j));
                } else {
                   Bn_posY(i,j,k) = dV(i,j,k)/(exp(dV(i,j,k)) - 1.0);
                   Bn_negY(i,j,k) = Bn_posY(i,j,k)*exp(dV(i,j,k));
                }
            }
        }
    }

}

void Continuity_n::Bernoulli_n_Z(const Eigen::Tensor<double, 3> &V_matrix)
{
    Eigen::Tensor<double, 3> dV(num_cell+1, num_cell+1, num_cell+1);
    dV.setZero();

    for (int i = 1; i < num_cell+1; i++)
       for (int j = 1; j < num_cell+1; j++)
           for (int k = 1; k < num_cell+1; k++)
                dV(i,j,k) =  V_matrix(i,j,k)-V_matrix(i,j,k-1);

    for (int i = 1; i < num_cell+1; i++) {
        for (int j = 1; j < num_cell+1; j++) {
            for (int k = 1; k < num_cell+1; k++) {
                if (abs(dV(i,j,k)) < 1e-13) {        //to prevent blowup due  to 0 denominator
                    Bn_posZ(i,j,k) = 1;//1 - dV(i,j)/2. + (dV(i,j)*dV(i,j))/12. - pow(dV(i,j), 4)/720.;
                    Bn_negZ(i,j,k) =  1;//Bn_posZ(i,j)*exp(dV(i,j));
                } else {
                   Bn_posZ(i,j,k) = dV(i,j,k)/(exp(dV(i,j,k)) - 1.0);
                   Bn_negZ(i,j,k) = Bn_posZ(i,j,k)*exp(dV(i,j,k));
                }
            }
        }
    }

}

//----------------------------------
//I THINK I DON'T NEED THIS ANYWHERE EXCEPT FOR WRITING TO FILE--> so move this conversion to utilities write to file function.
/*
void Continuity_n::to_matrix(const std::vector<double> &n)
{
    for (int index = 1; index <= num_elements; index++) {
        int i = index % N;    //this gives the i value for matrix
        if (i == 0) i = N;

        int j = 1 + static_cast<int>(floor((index-1)/N));  // j value for matrix

        n_matrix(i,j) = n[index];
    }
    for (int j = 1; j <= N; j++) {
        n_matrix(0, j) = n_leftBC[j];
        n_matrix(num_cell, j) = n_rightBC[j];
    }
    for (int i = 0; i <= num_cell; i++) {  //bottom BC's go all the way accross, including corners
        n_matrix(i, 0) = n_bottomBC[i];
        n_matrix(i, num_cell) = n_topBC[i];
    }

}
*/

void Continuity_n::calculate_currents()
{
    for (int i = 1; i < num_cell; i++) {
        for (int j = 1; j < num_cell; j++) {
            Jn_Z(i,j) =  J_coeff * n_mob(i,j) * (n_matrix(i,j)*Bn_posZ(i,j) - n_matrix(i,j-1)*Bn_negZ(i,j));
            Jn_X(i,j) =  J_coeff * n_mob(i,j) * (n_matrix(i,j)*Bn_posX(i,j) - n_matrix(i-1,j)*Bn_negX(i,j));
        }
    }
}
