#include "poisson.h"
#include <iostream>

Poisson::Poisson(const Parameters &params)
{
    main_diag.resize(num_elements+1);
    upper_diag.resize(num_elements);
    lower_diag.resize(num_elements);
    far_lower_diag.resize(num_elements-N+1);
    far_upper_diag.resize(num_elements-N+1);
    rhs.resize(num_cell);

    epsilon.resize(num_cell+1, num_cell+1);  //Eigen matrix object


    CV = params.N_dos*params.dx*params.dx*q/(epsilon_0*Vt);
    N = params.num_cell -1;  //for convenience define this --> is the number of points in 1D inside the device
    num_elements = params.num_elements;
    num_cell = params.num_cell;


    //FILL THE top and bottom BC's here


} //constructor)


void Poisson::setup_matrix()  //Note: this is on purpose different than the setup_eqn used for Continuity eqn's, b/c I need to setup matrix only once
{
    set_main_diag();
    set_lower_diag();
    set_upper_diag();
    set_far_lower_diag();
    set_far_upper_diag();
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

void Poisson::set_rhs(const std::vector<double> &n, const std::vector<double> &p, std::vector<double> &V_leftBC, std::vector<double> &V_rightBC, std::vector<double> &V_bottomBC, std::vector<double> &V_topBC)
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
                }else{
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
        }
