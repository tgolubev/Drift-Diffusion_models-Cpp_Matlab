#include "poisson.h"
#include <iostream>

Poisson::Poisson(const Parameters &params, const Eigen::MatrixXd &netcharge)
{
    CV = params.N_dos*params.dx*params.dx*q/(epsilon_0*Vt);
    N = params.num_cell -1;  //for convenience define this --> is the number of points in 1D inside the device
    num_elements = params.num_elements;
    num_cell = params.num_cell;
    V_matrix = Eigen::MatrixXd::Zero(num_cell+1, num_cell+1);    //useful for calculating currents at end of each Va

    main_diag.resize(num_elements+1);
    upper_diag.resize(num_elements);
    lower_diag.resize(num_elements);
    far_lower_diag.resize(num_elements-N+1);
    far_upper_diag.resize(num_elements-N+1);
    rhs.resize(num_elements+1);  //+1 b/c I am filling from index 1

    V_leftBC.resize(num_cell+1);
    V_rightBC.resize(num_cell+1);
    V_bottomBC.resize(num_cell+1);
    V_topBC.resize(num_cell+1);

    //MUST FILL WITH THE VALUES OF EPSILON!!  WILL NEED TO MODIFY THIS WHEN HAVE SPACE VARYING EPSILON
    epsilon =  params.eps_active*Eigen::MatrixXd::Ones(num_cell+1, num_cell+1);  //can fill with constant epsilon, like this!

    //allocate memory for the sparse matrix and rhs vector (Eig object)
    sp_matrix.resize(num_elements, num_elements);
    VecXd_rhs.resize(num_elements);   //only num_elements, b/c filling from index 0 (necessary for the sparse solver)

} //constructor)

//-------------------------------------------------------
//Set BC's functions


void Poisson::set_V_topBC(const Parameters &params, double Va)
{
    for (int i = 0; i <= num_cell; i++)
        V_topBC[i] = (params.Vbi-Va)/(2*Vt) - params.phi_c/Vt;
}

void Poisson::set_V_bottomBC(const Parameters &params, double Va)
{
    for (int i = 0; i <= num_cell; i++)
        V_bottomBC[i] = -((params.Vbi-Va)/(2*Vt) - params.phi_a/Vt);
}

void Poisson::set_V_leftBC(const std::vector<double> &V)
{
    for(int i = 1;i<=N;i++)
        V_leftBC[i] = V[(i-1)*N + 1];

}

void Poisson::set_V_rightBC(const std::vector<double> &V)
{
    for(int i = 1;i<=N;i++)
        V_rightBC[i] = V[i*N];
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


//---------------Setup AV diagonals (Poisson solve)---------------------------------------------------------------

//far lower diagonal
void Poisson::set_far_lower_diag(){
    for(int index = 1; index<= N*(N-1); index++){
        int i = index % N;
        if(i==0) i=N;
        int j = 1 + static_cast<int>(floor((index-1)/N));

        far_lower_diag[index] = -(epsilon(i,j) + epsilon(i+1,j))/2.;
    }
}


//lower diag
void Poisson::set_lower_diag(){
    //int num_elements = lower_diag.size() -1;

    for (int index = 1; index<=num_elements-1;index++){  //      %this is the lower diagonal (below main diagonal) (1st element corresponds to 2nd row)
        int i = index % N;        // %this is x index of V which element corresponds to (note if this = 0, means these are the elements which are 0);
        int j = 1 + static_cast<int>(floor((index-1)/N));

        if(index % N == 0)
            lower_diag[index] = 0; //  %these are the elements at subblock corners
        else
            lower_diag[index] = -(epsilon(i,j) + epsilon(i,j+1))/2.;
    }

}


//main diagonal
void Poisson::set_main_diag(){
    //int num_elements = main_diag.size() - 1;

    for (int index =  1; index<=num_elements;index++){//      %main diagonal
        int i = index % N;
        if(i ==0)        //        %the multiples of N correspond to last index
            i = N;
        int j = 1 + static_cast<int>(floor((index-1)/N));

        main_diag[index] = (epsilon(i+1,j) + epsilon(i+1,j+1))/2. + (epsilon(i,j) + epsilon(i,j+1))/2. + (epsilon(i,j+1) + epsilon(i+1,j+1))/2. + (epsilon(i,j) + epsilon(i+1,j))/2.;
    }
}

//upper diagonal
void Poisson::set_upper_diag(){
    //int num_elements = upper_diag.size() -1;  //-1 since vectors fill from 0, but I'm indexing from 1
    //NOTE: num_elements here, is actually # of elements in the off-diagonal

    for (int index = 1; index <=num_elements-1;index++) {  //      %main uppper diagonal, matlab fills this from the bottom (so i = 2 corresponds to 1st row in matrix)
        int i = index % N;
        int j = 1 + static_cast<int>(floor((index-1)/N));

        if(index  % N ==0)
            upper_diag[index] = 0;
        else
            upper_diag[index] =  -(epsilon(i+1,j) + epsilon(i+1,j+1))/2.;
   }
}


//far upper diagonal
void Poisson::set_far_upper_diag(){

    for (int index = 1; index <= num_elements-N; index++) { //      %far upper diagonal, matlab fills from bottom, so this starts at 1+N (since 1st element is in the 2nd subblock of matrix)
        int i = index % N;
        if(i ==0)      //          %the multiples of N correspond to last index
            i = N;
        int j = 1 + static_cast<int>(floor((index-N)/N));

         far_upper_diag[index] = -(epsilon(i,j+1) + epsilon(i+1,j+1))/2.;         //    %1st element corresponds to 1st row.   this has N^2 -N elements
    }
}

void Poisson::set_rhs(const Eigen::MatrixXd &netcharge)
{
    //setup rhs of Poisson eqn.
    int index2 = 0;
    for(int j = 1;j<=N;j++){
        if(j==1){
            for(int i = 1;i<=N;i++){
                index2++;                //THIS COULD BE MODIFIES TO a switch ,case statement--> might be cleaner
                if(i==1){
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*(V_leftBC[1] + V_bottomBC[i]);
                }else if(i == N)
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*(V_rightBC[1] + V_bottomBC[i]);
                else
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*V_bottomBC[i];
            }
        }else if(j==N){
            for(int i = 1; i<=N;i++){
                index2++;
                if(i==1)
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*(V_leftBC[N] + V_topBC[i]);
                else if(i == N)
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*(V_rightBC[N] + V_topBC[i]);
                else
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*V_topBC[i];
            }
        }else{  //these seems ok
            for(int i = 1;i<=N;i++){
                index2++;

                if(i==1)
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*V_leftBC[j];
                else if(i == N)
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*V_rightBC[j];
                else
                    rhs[index2] = netcharge(i,j);

            }
        }
    }


    //set up VectorXd Eigen vector object for sparse solver
    for (int i = 1; i<=num_elements; i++) {
        VecXd_rhs(i-1) = rhs[i];   //fill VectorXd  rhs of the equation
    }

}

//-----------------------------------------
void Poisson::to_matrix(const std::vector<double> &V)
{
    for (int index = 1; index <= num_elements; index++) {
        int i = index % N;    //this gives the i value for matrix
        if (i == 0) i = N;

        int j = 1 + static_cast<int>(floor((index-1)/N));  // j value for matrix

        V_matrix(i,j) = V[index];
    }

    for (int i = 0; i <=num_cell; i++) {
        V_matrix(i, 0) = V_bottomBC[i];
        V_matrix(i,num_cell) = V_topBC[i];
    }

    for (int j = 1; j < num_cell; j++) {   //don't need to set j = 0 and j = num_cell elements, b/c already set when apply top and bottom BC's (are corners)
        V_matrix(0, j) = V_leftBC[j];
        V_matrix(num_cell, j) = V_rightBC[j];
    }

}
