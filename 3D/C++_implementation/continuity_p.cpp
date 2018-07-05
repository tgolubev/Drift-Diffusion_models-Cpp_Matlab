#include "continuity_p.h"

Continuity_p::Continuity_p(const Parameters &params)
{
    num_elements = params.num_elements;
    N = params.num_cell - 1;
    num_cell = params.num_cell;
    p_matrix = Eigen::MatrixXd::Zero(num_cell+1, num_cell+1);

    main_diag.resize(num_elements+1);
    upper_diag.resize(num_elements+1);
    lower_diag.resize(num_elements+1);
    far_lower_diag.resize(num_elements+1);
    far_upper_diag.resize(num_elements+1);
    rhs.resize(num_elements+1);  //+1 b/c I am filling from index 1

    Bp_posX = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
    Bp_negX = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
    Bp_posY = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
    Bp_negY = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
    Bp_posZ = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
    Bp_negZ = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);

    values = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
    values2 = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
    values3 = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);
    values4 = Eigen::Tensor<double, 3> (num_cell+1, num_cell+1, num_cell+1);

    Jp_Z.resize(num_cell+1, num_cell+1);
    Jp_X.resize(num_cell+1, num_cell+1);

    p_bottomBC.resize(num_cell+1, num_cell+1);
    p_topBC.resize(num_cell+1, num_cell+1);
    p_leftBC_X.resize(num_cell+1, num_cell+1);
    p_rightBC_X.resize(num_cell+1, num_cell+1);
    p_leftBC_Y.resize(num_cell+1, num_cell+1);
    p_rightBC_Y.resize(num_cell+1, num_cell+1);

    J_coeff = (q*Vt*params.N_dos*params.mobil)/params.dx;

    //------------------------------------------------------------------------------------------
    //MUST FILL WITH THE VALUES OF p_mob!!  WILL NEED TO MODIFY THIS WHEN HAVE SPACE VARYING
    p_mob = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
    p_mob_avg_X = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
    p_mob_avg_Y = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
    p_mob_avg_Z = Eigen::Tensor<double, 3> (num_cell+2, num_cell+2, num_cell+2);
    p_mob.setConstant(params.p_mob_active/params.mobil);

    //Compute averaged mobilities
    for (int k = 1; k <= num_cell+1; k++) {
        for (int j = 1; j <= num_cell+1; j++) {
            for (int i = 1; i <= num_cell+1; i++) {
                p_mob_avg_X(i,j,k) = (p_mob(i,j,k) + p_mob(i,j+1,k) + p_mob(i,j,k+1) + p_mob(i,j+1,k+1))/4.;
                p_mob_avg_Y(i,j,k) = (p_mob(i,j,k) + p_mob(i+1,j,k) + p_mob(i,j,k+1) + p_mob(i+1,j,k+1))/4.;
                p_mob_avg_Z(i,j,k) = (p_mob(i,j,k) + p_mob(i+1,j,k) + p_mob(i,j+1,k) + p_mob(i+1,j+1,k))/4.;
            }
        }
    }

     //------------------------------------------------------------------------------------------
    Cp = (params.dx*params.dx)/(Vt*params.N_dos*params.mobil);  //can't use static, b/c dx wasn't defined as const, so at each initialization of Continuity_p object, new const will be made.

    //these BC's for now stay constant throughout simulation, so fill them once, upon Continuity_p object construction
    for (int j =  0; j <= num_cell; j++) {
        for (int i = 1; i <= num_cell+1; i++) {
            p_bottomBC(i,j) = params.N_HOMO*exp(-params.phi_a/Vt)/params.N_dos;
            p_topBC(i,j) = params.N_HOMO*exp(-(params.E_gap-params.phi_c)/Vt)/params.N_dos;
        }
    }

    //allocate memory for the sparse matrix and rhs vector (Eig object)
    sp_matrix.resize(num_elements, num_elements);
    VecXd_rhs.resize(num_elements);   //only num_elements, b/c filling from index 0 (necessary for the sparse solver)

    //setup the triplet list for sparse matrix
     triplet_list.resize(5*num_elements);   //approximate the size that need
}

//------------------------------------------------------------------
//Set BC's
void Continuity_p::set_p_leftBC_X(const std::vector<double> &p)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            p_leftBC_X(j,k) = p[index + (j-1)*N + 1];
        }
        index = index+N*N;  //brings us to next vertical subblock set
    }
}

void Continuity_p::set_p_rightBC_X(const std::vector<double> &p)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            p_rightBC_X(j,k) = p[index + j*N];
        }
        index = index+N*N;  //brings us to next vertical subblock set
    }
}

void Continuity_p::set_p_leftBC_Y(const std::vector<double> &p)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int i = 1; i <= N; i++) {
            p_leftBC_Y(i, k) = p[index + i];
        }
        index = index + N*N;
    }
}

void Continuity_p::set_p_rightBC_Y(const std::vector<double> &p)
{
    int index = 0;
    for (int k = 1; k <= N; k++) {
        for (int i = 1; i <= N; i++) {
            p_rightBC_Y(i, k) = p[index + i + N*N - N];
        }
        index = index + N*N;
    }
}


//Calculates Bernoulli fnc values, then sets the diagonals and rhs
void Continuity_p::setup_eqn(const Eigen::Tensor<double, 3> &V_matrix, const std::vector<double> &Up, const std::vector<double> &p)
{
    trp_cnt = 0;  //reset triplet count
    Bernoulli_p_X(V_matrix);
    Bernoulli_p_Y(V_matrix);
    Bernoulli_p_Z(V_matrix);

    set_far_lower_diag();
    set_lower_diag();
    set_main_lower_diag();
    set_main_diag();
    set_main_upper_diag();
    set_upper_diag();
    set_far_upper_diag();

    set_p_leftBC_X(p);
    set_p_rightBC_X(p);
    set_p_leftBC_Y(p);
    set_p_rightBC_Y(p);

    set_rhs(Up);

    sp_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());    //sp_matrix is our sparse matrix
}

//------------------------------Setup Ap diagonals----------------------------------------------------------------
void Continuity_p::set_far_lower_diag()
{
    values = -p_mob_avg_Z*Bp_posZ;  //element wise multiplication

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

void Continuity_p::set_lower_diag()
{
    values = -p_mob_avg_Y*Bp_posY;

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
void Continuity_p::set_main_lower_diag()
{
    values = -p_mob_avg_X*Bp_posX;

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


void Continuity_p::set_main_diag()
{
    values  = p_mob_avg_Z*Bp_negZ + p_mob_avg_Y*Bp_negY + p_mob_avg_X*Bp_negX;
    values2 = p_mob_avg_X*Bp_posX;
    values3 = p_mob_avg_Y*Bp_posY;
    values4 = p_mob_avg_Z*Bp_posZ;

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


void Continuity_p::set_main_upper_diag()
{
    values = -p_mob_avg_X*Bp_negX;

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


void Continuity_p::set_upper_diag()
{
    values = -p_mob_avg_Y*Bp_negY;

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


void Continuity_p::set_far_upper_diag()
{
  values = -p_mob_avg_Z*Bp_negZ;

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

//---------------------------

void Continuity_p::Bernoulli_p_X(const Eigen::Tensor<double, 3> &V_matrix)
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
                    Bp_posX(i,j,k) = 1;//1 - dV(i,j)/2. + (dV(i,j)*dV(i,j))/12. - pow(dV(i,j), 4)/720.;
                    Bp_negX(i,j,k) =  1;//Bp_posX(i,j)*exp(dV(i,j));
                } else {
                    Bp_posX(i,j,k) = dV(i,j,k)/(exp(dV(i,j,k)) - 1.0);
                    Bp_negX(i,j,k) = Bp_posX(i,j,k)*exp(dV(i,j,k));
                }
            }
        }
    }
}

void Continuity_p::Bernoulli_p_Y(const Eigen::Tensor<double, 3> &V_matrix)
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
                    Bp_posY(i,j,k) = 1;//1 - dV(i,j)/2. + (dV(i,j)*dV(i,j))/12. - pow(dV(i,j), 4)/720.;
                    Bp_negY(i,j,k) =  1;//Bp_posZ(i,j)*exp(dV(i,j));
                } else {
                   Bp_posY(i,j,k) = dV(i,j,k)/(exp(dV(i,j,k)) - 1.0);
                   Bp_negY(i,j,k) = Bp_posY(i,j,k)*exp(dV(i,j,k));
                }
            }
        }
    }

}

void Continuity_p::Bernoulli_p_Z(const Eigen::Tensor<double, 3> &V_matrix)
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
                    Bp_posZ(i,j,k) = 1;//1 - dV(i,j)/2. + (dV(i,j)*dV(i,j))/12. - pow(dV(i,j), 4)/720.;
                    Bp_negZ(i,j,k) =  1;//Bp_posZ(i,j)*exp(dV(i,j));
                } else {
                   Bp_posZ(i,j,k) = dV(i,j,k)/(exp(dV(i,j,k)) - 1.0);
                   Bp_negZ(i,j,k) = Bp_posZ(i,j,k)*exp(dV(i,j,k));
                }
            }
        }
    }

}


//--------------------------------------------------------------------------------
void Continuity_p::set_rhs(const std::vector<double> &Up)
{
    //calculate main part here
    for (int i = 1; i <= num_elements; i++)
        rhs[i] = Cp*Up[i];

    //add on BC's
    int index = 0;
    for (int k = 1; k <= N; k++) {
        if (k == 1) {
            for (int j = 1; j <= N; j++) {
                if (j == 1)  {//different for 1st subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)      //1st element has 2 BC's
                            rhs[index] += p_mob(i,j,k)*p_leftBC_X(1,1)*Bp_posX(i,j,k) + p_mob(i,j,k)*p_leftBC_Y(1,1)*Bp_posY(i,j,k) + p_mob(i,j,k)*p_bottomBC(i,j)*Bp_posZ(i,j,k);
                        else if (i == N)
                            rhs[index] += p_mob(N+1, 1, 1)*p_rightBC_X(1,1)*Bp_negX(i+1,j,k) + p_mob(N, 0, 1)*p_leftBC_Y(N,1)*Bp_posY(i,j,k) + p_mob(N, 1, 0)*p_bottomBC(i,j)*Bp_posZ(i,j,k);
                        else //middle elements in 1st subblock
                            rhs[index] += p_mob(i, 0, 1)*p_leftBC_Y(i, 1)*Bp_posY(i,j,k) + p_mob(i, 1, 0)*p_bottomBC(i,j)*Bp_posZ(i,j,k);
                    }
                } else if (j == N) {      //different for last subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)  //1st element has 2 BC's
                            rhs[index] += p_mob(0, N, 1)*p_leftBC_X(N,1)*Bp_posX(i,j,k) + p_mob(1, N+1, 1)*p_rightBC_Y(1,1)*Bp_negY(i,j+1,k) + p_mob(1, N, 0)*p_bottomBC(i,j)*Bp_posZ(i,j,k);
                        else if (i==N)
                            rhs[index] += p_mob(N+1, N, 1)*p_rightBC_X(N,1)*Bp_negX(i+1,j,k) + p_mob(N, N+1, 1)*p_rightBC_Y(N,1)*Bp_negY(i,j+1,k) + p_mob(N,N, 0)*p_bottomBC(i,j)*Bp_posZ(i,j,k);
                        else //inner rows of Nth subblock
                            rhs[index] += p_mob(i,N+1, 1)*p_rightBC_Y(i,1)*Bp_negY(i,j+1,k) + p_mob(i, N, 0)*p_bottomBC(i,j)*Bp_posZ(i,j,k);
                    }
                } else {     //interior subblocks of 1st (k = 1) subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)
                            rhs[index] += p_mob(0, j, 1)*p_leftBC_X(j,1)*Bp_posX(i,j,k) + p_mob(1, j, 0)*p_bottomBC(i,j)*Bp_posZ(i,j,k);
                        else if (i == N)
                            rhs[index] += p_mob(N+1, j, 1)*p_rightBC_X(j,1)*Bp_negX(i+1,j,k) + p_mob(N, j, 0)*p_bottomBC(i,j)*Bp_posZ(i,j,k);
                        else
                            rhs[index] += p_mob(i, j, 0)*p_bottomBC(i,j)*Bp_posZ(i,j,k);
                    }
                }
            }
        } else if (k == N) {
            for (int j = 1; j <= N; j++) {
                if (j == 1)  {//different for 1st subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)      //1st element has 2 BC's
                            rhs[index] += p_mob(0,1,N)*p_leftBC_X(1,N)*Bp_posX(i,j,k) + p_mob(1,0,N)*p_leftBC_Y(1,N)*Bp_posY(i,j,k) + p_mob(1,1,N+1)*p_topBC(i,j)*Bp_negZ(i,j,k+1);
                        else if (i == N)
                            rhs[index] += p_mob(N+1, 1, N)*p_rightBC_X(1,N)*Bp_negX(i+1,j,k) + p_mob(N, 0, N)*p_leftBC_Y(N,N)*Bp_posY(i,j,k) + p_mob(N, 1, N+1)*p_topBC(i,j)*Bp_negZ(i,j,k+1);
                        else //middle elements in 1st subblock
                            rhs[index] += p_mob(i, 0, N)*p_leftBC_Y(i, N)*Bp_posY(i,j,k) + p_mob(i, 1, N+1)*p_topBC(i,j)*Bp_negZ(i,j,k+1);
                    }
                } else if (j == N) {      //different for last subblock within k=N subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)  //1st element has 2 BC's
                            rhs[index] += p_mob(0,N, N)*p_leftBC_X(N,N)*Bp_posX(i,j,k) + p_mob(1, N+1, N)*p_rightBC_Y(1,N)*Bp_negY(i,j+1,k) + p_mob(1, N, N+1)*p_topBC(i,j)*Bp_negZ(i,j,k+1);
                        else if (i==N)
                            rhs[index] += p_mob(N+1,N, N)*p_rightBC_X(N,N)*Bp_negX(i+1,j,k) + p_mob(N, N+1, N)*p_rightBC_Y(N,N)*Bp_negY(i,j+1,k) + p_mob(N,N, N+1)*p_topBC(i,j)*Bp_negZ(i,j,k+1);
                        else //inner rows of Nth subblock
                            rhs[index] += p_mob(i,N+1, N)*p_rightBC_Y(i,N)*Bp_negY(i,j+1,k) + p_mob(i, N, N+1)*p_topBC(i,j)*Bp_negZ(i,j,k+1);
                    }
                } else {    //interior subblocks of last (k = N) subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)
                            rhs[index] += p_mob(0, j, N)*p_leftBC_X(j,N)*Bp_posX(i,j,k) + p_mob(1, j, N+1)*p_topBC(i,j)*Bp_negZ(i,j,k+1);
                        else if (i == N)
                            rhs[index] += p_mob(N+1, j, N)*p_rightBC_X(j,N)*Bp_negX(i+1,j,k) + p_mob(N, j, N+1)*p_topBC(i,j)*Bp_negZ(i,j,k+1);
                        else
                            rhs[index] += p_mob(i, j, N+1)*p_topBC(i,j)*Bp_negZ(i,j,k+1);
                    }
                }
            }
        } else {  //interior subblock groups (k=2:N-1)
            for (int j = 1; j <= N; j++) {
                if (j == 1)  {//different for 1st subblock
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)      //1st element has 2 BC's
                            rhs[index] += p_mob(0, 1, k)*p_leftBC_X(1,k)*Bp_posX(i,j,k)  + p_mob(1, 0, k)*p_leftBC_Y(1,k)*Bp_posY(i,j,k);
                        else if (i == N)
                            rhs[index] += p_mob(N+1, 1, k)*p_rightBC_X(1,k)*Bp_negX(i+1,j,k) + p_mob(N, 0, k)*p_leftBC_Y(N,k)*Bp_posY(i,j,k);
                        else //middle elements in 1st subblock
                            rhs[index] += p_mob(i, 0, k)*p_leftBC_Y(i, k)*Bp_posY(i,j,k);
                    }
                } else if (j == N) {      //different for last subblock within k=N subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)  //1st element has 2 BC's
                            rhs[index] += p_mob(0,N, k)*p_leftBC_X(N,k)*Bp_posX(i,j,k) + p_mob(1, N+1, k)*p_rightBC_Y(1,k)*Bp_negY(i,j+1,k);
                        else if (i==N)
                            rhs[index] += p_mob(N+1,N, k)*p_rightBC_X(N,k)*Bp_negX(i+1,j,k) + p_mob(N, N+1, k)*p_rightBC_Y(N,k)*Bp_negY(i,j+1,k);
                        else //inner rows of Nth subblock
                            rhs[index] += p_mob(i,N+1, k)*p_rightBC_Y(i,k)*Bp_negY(i,j+1,k);
                    }
                } else {    //interior subblocks of last (k = N) subblock group
                    for (int i = 1; i <= N; i++) {
                        index++;
                        if (i == 1)
                            rhs[index] += p_mob(0, j, k)*p_leftBC_X(j,k)*Bp_posX(i,j,k) ;
                        else if (i == N)
                            rhs[index] += p_mob(N+1, j, k)*p_rightBC_X(j,k)*Bp_negX(i+1,j,k);
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

//----------------------------------------------
//I THINK I DON'T NEED THIS ANYWHERE EXCEPT FOR WRITING TO FILE--> so move this conversion to utilities write to file function.
/*
void Continuity_p::to_matrix(const std::vector<double> &p)
{
    for (int index = 1; index <= num_elements; index++) {
        int i = index % N;    //this gives the i value for matrix
        if (i == 0) i = N;

        int j = 1 + static_cast<int>(floor((index-1)/N));  // j value for matrix

        p_matrix(i,j) = p[index];
    }
    for (int j = 1; j <= N; j++) {
        p_matrix(0, j) = p_leftBC[j];
        p_matrix(num_cell, j) = p_rightBC[j];
    }
    for (int i = 0; i <= num_cell; i++) {  //bottom BC's go all the way accross, including corners
        p_matrix(i, 0) = p_bottomBC[i];
        p_matrix(i, num_cell) = p_topBC[i];
    }

}
*/


void Continuity_p::calculate_currents()
{
    for (int i = 1; i <= num_cell; i++) {
        for (int j = 1; j < num_cell; j++) {
            Jp_Z(i,j) = -J_coeff * p_mob(i,j) * (p_matrix(i,j)*Bp_negZ(i,j) - p_matrix(i,j-1)*Bp_posZ(i,j));
            Jp_X(i,j) = -J_coeff * p_mob(i,j) * (p_matrix(i,j)*Bp_negX(i,j) - p_matrix(i-1,j)*Bp_posX(i,j));
        }
    }

}
