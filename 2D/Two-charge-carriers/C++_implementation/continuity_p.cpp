#include "continuity_p.h"

Continuity_p::Continuity_p(const Parameters &params)
{
    num_elements = params.num_elements;
    N = params.num_cell - 1;
    num_cell = params.num_cell;

    main_diag.resize(num_elements+1);
    upper_diag.resize(num_elements);
    lower_diag.resize(num_elements);
    far_lower_diag.resize(num_elements-N+1);
    far_upper_diag.resize(num_elements-N+1);
    rhs.resize(num_elements+1);  //+1 b/c I am filling from index 1

    Bp_posX.resize(num_cell+1, num_cell+1);  //allocate memory for the matrix object
    Bp_negX.resize(num_cell+1, num_cell+1);
    Bp_posZ.resize(num_cell+1, num_cell+1);
    Bp_negZ.resize(num_cell+1, num_cell+1);

    p_bottomBC.resize(N+1);
    p_topBC.resize(N+1);
    p_leftBC.resize(N+1);
    p_rightBC.resize(N+1);

    // //MUST FILL WITH THE VALUES OF p_mob!!  WILL NEED TO MODIFY THIS WHEN HAVE SPACE VARYING
    p_mob = (params.p_mob_active/params.mobil)*Eigen::MatrixXd::Ones(num_cell+1, num_cell+1);
    //p_mob.resize(num_cell+1, num_cell+1);  //note: p_mob is an Eigen Matrix object...

    Cp = params.dx*params.dx/(Vt*params.N_dos*params.mobil);  //can't use static, b/c dx wasn't defined as const, so at each initialization of Continuity_p object, new const will be made.

    //these BC's for now stay constant throughout simulation, so fill them once, upon Continuity_n object construction

    for (int j =  0; j <= N; j++) {
        p_bottomBC[j] = params.N_HOMO*exp(-params.phi_a/Vt)/params.N_dos;
        p_topBC[j] = params.N_HOMO*exp(-(params.E_gap-params.phi_c)/Vt)/params.N_dos;
    }

    //allocate memory for the sparse matrix and rhs vector (Eig object)
    sp_matrix.resize(num_elements, num_elements);
    VecXd_rhs.resize(num_elements);   //only num_elements, b/c filling from index 0 (necessary for the sparse solver)

}

//------------------------------------------------------------------
//Set BC's
void Continuity_p::set_p_leftBC(const std::vector<double> &p)
{
    for (int j = 1; j <= N; j++) {
        p_leftBC[j] = p[(j-1)*N + 1];
    }
}

void Continuity_p::set_p_rightBC(const std::vector<double> &p)
{
    for (int j = 1; j <= N; j++) {
         p_rightBC[j]= p[j*N];
    }
}




//Calculates Bernoulli fnc values, then sets the diagonals and rhs
void Continuity_p::setup_eqn(const Eigen::MatrixXd &V_matrix, const Eigen::MatrixXd &Up_matrix, const std::vector<double> &p)
{
    Bernoulli_p_X(V_matrix);
    Bernoulli_p_Z(V_matrix);
    set_far_lower_diag();
    set_lower_diag();
    set_main_diag();
    set_upper_diag();
    set_far_upper_diag();
    set_p_leftBC(p);
    set_p_rightBC(p);
    set_rhs(Up_matrix);

    typedef Eigen::Triplet<double> Trp;

    //THIS COULD BE MOVED TO UTILITIES..., B/C IS BASICALLY THE SAME for the 2 continuity and poisson equations.

    //generate triplets for Eigen sparse matrix
    //setup the triplet list for sparse matrix
     std::vector<Trp> triplet_list(5*num_elements);   //approximate the size that need         // list of non-zeros coefficients in triplet form(row index, column index, value)
     int trp_cnt = 0;
     for(int i = 1; i<= num_elements; i++){
           triplet_list[trp_cnt] = {i-1, i-1, main_diag[i]};
           trp_cnt++;
         //triplet_list.push_back(Trp(i-1,i-1,main_diag[i]));;   //triplets for main diag
     }
     for(int i = 1;i< upper_diag.size();i++){
         triplet_list[trp_cnt] = {i-1, i, upper_diag[i]};
         trp_cnt++;
         //triplet_list.push_back(Trp(i-1, i, upper_diag[i]));  //triplets for upper diagonal
      }
      for(int i = 1;i< lower_diag.size();i++){
          triplet_list[trp_cnt] = {i, i-1, lower_diag[i]};
          trp_cnt++;
         // triplet_list.push_back(Trp(i, i-1, lower_diag[i]));
      }
      for(int i = 1;i< far_upper_diag.size();i++){
          triplet_list[trp_cnt] = {i-1, i-1+N, far_upper_diag[i]};
          trp_cnt++;
          //triplet_list.push_back(Trp(i-1, i-1+N, far_upper_diag[i]));
          triplet_list[trp_cnt] = {i-1+N, i-1, far_lower_diag[i]};
          trp_cnt++;
          //triplet_list.push_back(Trp(i-1+N, i-1, far_lower_diag[i]));
       }

     sp_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());    //sp_matrix is our sparse matrix
}

//------------------------------Setup Ap diagonals----------------------------------------------------------------
void Continuity_p::set_far_lower_diag()
{
    //Lowest diagonal: corresponds to V(i, j-1)
    for (int index = 1; index <=N*(N-1); index++) {      //(1st element corresponds to Nth row  (number of elements = N*(N-1)
        int i = index % N;

        if (i ==0)                //the multiples of N correspond to last index
            i = N;

        int j = 1 + static_cast<int>(floor((index-1)/N));    //this is the y index of V which element corresponds to. 1+ floor(index/4)determines which subblock this corresponds to and thus determines j, since the j's for each subblock are all the same.

        far_lower_diag[index] = -((p_mob(i,j) + p_mob(i+1, j))/2.)*Bp_posZ(i,j);
    }
}


//main lower diag
void Continuity_p::set_lower_diag()
{
    for (int index = 1; index <= num_elements-1; index++) {
        int i = index % N;
        int j = 1 + static_cast<int>(floor((index-1)/N));

        if (index % N == 0)
            lower_diag[index] = 0;      //probably don't need explicitely fill here, since auto intialized to 0
        else
            lower_diag[index] = -((p_mob(i,j) + p_mob(i,j+1))/2.)*Bp_posX(i,j);
    }
}


void Continuity_p::set_main_diag()
{
    for (int index = 1; index <= num_elements; index++) {
        int i = index % N;

        if (i == 0)
            i = N;

        int j = 1 + static_cast<int>(floor((index-1)/N));

        main_diag[index] = ((p_mob(i,j) + p_mob(i,j+1))/2.)*Bp_negX(i,j) + ((p_mob(i+1,j) + p_mob(i+1,j+1))/2.)*Bp_posX(i+1,j) + ((p_mob(i,j) + p_mob(i+1,j))/2.)*Bp_negZ(i,j) + ((p_mob(i,j+1) + p_mob(i+1,j+1))/2.)*Bp_posZ(i,j+1);
    }
}


void Continuity_p::set_upper_diag()
{
    for (int index = 1; index <= num_elements-1; index++) {
        int i = index % N;
        int j = 1 + static_cast<int>(floor((index-1)/N));

        if ((index-1 % N) == 0)
            upper_diag[index] = 0;
        else
            upper_diag[index] = -((p_mob(i+1,j) + p_mob(i+1,j+1))/2.)*Bp_negX(i+1,j);
    }
}


void Continuity_p::set_far_upper_diag()
{
    for (int index = 1; index <= num_elements-N; index++) {
        int i = index % N;

        if(i == 0)
            i = N;

        int j = 1 + static_cast<int>(floor((index-N)/N));

        far_upper_diag[index] = -((p_mob(i,j+1) + p_mob(i+1,j+1))/2.)*Bp_negZ(i,j+1);
    }
}


//---------------------------

void Continuity_p::Bernoulli_p_X(const Eigen::MatrixXd &V_matrix)
{
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(num_cell+1,num_cell+1);

    for (int i = 1; i < num_cell+1; i++)
       for (int j = 1; j < num_cell+1; j++)
           dV(i,j) =  V_matrix(i,j)-V_matrix(i-1,j);

    for (int i = 1; i < num_cell+1; i++) {
        for (int j = 1; j < num_cell+1; j++) {
            if (abs(dV(i,j)) < 1e-13) {        //to prevent blowup due  to 0 denominator
                Bp_posX(i,j) = 1 - dV(i,j)/2. + (dV(i,j)*dV(i,j))/12. - pow(dV(i,j), 4)/720.;
                Bp_negX(i,j) =  Bp_posX(i,j)*exp(dV(i));
            } else {
               Bp_posX(i,j) = dV(i,j)/(exp(dV(i,j)) - 1.0);
               Bp_negX(i,j) = Bp_posX(i,j)*exp(dV(i,j));
            }
        }
    }
}

void Continuity_p::Bernoulli_p_Z(const Eigen::MatrixXd &V_matrix)
{
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(num_cell+1,num_cell+1);

    for (int i = 1; i < num_cell+1; i++)
       for (int j = 1; j < num_cell+1; j++)
           dV(i,j) =  V_matrix(i,j)-V_matrix(i,j-1);

    for (int i = 1; i < num_cell+1; i++) {
        for (int j = 1; j < num_cell+1; j++) {
            if (abs(dV(i,j)) < 1e-13) {        //to prevent blowup due  to 0 denominator
                Bp_posZ(i,j) = 1 - dV(i,j)/2. + (dV(i,j)*dV(i,j))/12. - pow(dV(i,j), 4)/720.;
                Bp_negZ(i,j) =  Bp_posZ(i,j)*exp(dV(i));
            } else {
               Bp_posZ(i,j) = dV(i,j)/(exp(dV(i,j)) - 1.0);
               Bp_negZ(i,j) = Bp_posZ(i,j)*exp(dV(i,j));
            }
        }
    }
}

//--------------------------------------------------------------------------------
void Continuity_p::set_rhs(const Eigen::MatrixXd &Up_matrix)
{
    int index = 0;

    for (int j = 1; j <= N; j++) {
        if (j ==1)  {//different for 1st subblock
            for (int i = 1; i <= N; i++) {
                index++;
                if (i==1)     //1st element has 2 BC's
                    rhs[index] = Cp*Up_matrix(i,j) + p_mob(i,j)*(Bp_posX(i,j)*p_leftBC[1] + Bp_posZ(i,j)*p_bottomBC[i]);  //NOTE: rhs is +Cp*Up_matrix, b/c diagonal elements are + here, flipped sign from 1D version
                else if (i==N)
                    rhs[index] = Cp*Up_matrix(i,j) + p_mob(i,j)*(Bp_posZ(i,j)*p_bottomBC[i] + Bp_negX(i+1,j)*p_rightBC[1]);
                else
                    rhs[index] = Cp*Up_matrix(i,j) + p_mob(i,j)*Bp_posZ(i,j)*p_bottomBC[i];
            }
        } else if (j == N) {      //different for last subblock
            for (int i = 1; i <= N; i++) {
                index++;
                if (i==1)  //1st element has 2 BC's
                    rhs[index] = Cp*Up_matrix(i,j) + p_mob(i,j)*(Bp_posX(i,j)*p_leftBC[N] + Bp_negZ(i,j+1)*p_topBC[i]);
                else if (i==N)
                        rhs[index] = Cp*Up_matrix(i,j) + p_mob(i,j)*(Bp_negX(i+1,j)*p_rightBC[N] + Bp_negZ(i,j+1)*p_topBC[i]);
                else
                rhs[index] = Cp*Up_matrix(i,j) + p_mob(i,j)*Bp_negZ(i,j+1)*p_topBC[i];
            }
        } else {     //interior subblocks
            for (int i = 1; i <= N; i++) {
                index++;
                if(i==1)
                    rhs[index] = Cp*Up_matrix(i,j) + p_mob(i,j)*Bp_posX(i,j)*p_leftBC[j];
                else if(i==N)
                        rhs[index] = Cp*Up_matrix(i,j) + p_mob(i,j)*Bp_negX(i+1,j)*p_rightBC[j];
                else
                rhs[index] = Cp*Up_matrix(i,j);
            }
        }
    }

    //set up VectorXd Eigen vector object for sparse solver
    for (int i = 1; i<=num_elements; i++) {
        VecXd_rhs(i-1) = rhs[i];   //fill VectorXd  rhs of the equation
    }
}
