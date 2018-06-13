#include "continuity_n.h"

Continuity_n::Continuity_n(const Parameters &params)
{
   num_cell = params.num_cell;
   main_diag.resize(num_cell);
   upper_diag.resize(num_cell-1);
   lower_diag.resize(num_cell-1);
   rhs.resize(num_cell);

   Bn_posX.resize(num_cell+1, num_cell+1);  //allocate memory for the matrix object
   Bn_negX.resize(num_cell+1, num_cell+1);
   Bn_posZ.resize(num_cell+1, num_cell+1);  //allocate memory for the matrix object
   Bn_negZ.resize(num_cell+1, num_cell+1);

   n_mob.resize(num_cell+1, num_cell+1);  //note: n_mob is an Eigen Matrix object...

   Cn = params.dx*params.dx/(Vt*params.N_dos*params.mobil);

   for (int j =  0; j <= N; j++) {
      n_bottomBC[j] = params.N_LUMO*exp(-(params.E_gap-params.phi_a)/Vt)/params.N_dos;
      n_topBC[j] = params.N_LUMO*exp(-params.phi_c/Vt)/params.N_dos;
   }

   num_elements = params.num_elements;
   N = params.num_cell - 1;

}

//Calculates Bernoulli fnc values, then sets the diagonals and rhs
//use the V_matrix for setup, to be able to write equations in terms of (x,z) coordingates
void Continuity_n::setup_eqn(const Eigen::MatrixXd &V_matrix, const Eigen::MatrixXd &Un_matrix)
{
    Bernoulli_n_X(V_matrix);
    Bernoulli_n_Z(V_matrix);
    set_main_diag();
    set_upper_diag();
    set_lower_diag();
    set_rhs(Un_matrix);

}

//-------------------------------Setup An diagonals (Continuity/drift-diffusion solve)-----------------------------

void Continuity_n::set_far_lower_diag()
{
    //Lowest diagonal: corresponds to V(i, j-1)
    for (int index = 1; index <=N*(N-1); index++) {      //(1st element corresponds to Nth row  (number of elements = N*(N-1)
        int i = index % N;

        if (i ==0)                //the multiples of N correspond to last index
            i = N;

        int j = 1 + floor((index-1)/N);    //this is the y index of V which element corresponds to. 1+ floor(index/4)determines which subblock this corresponds to and thus determines j, since the j's for each subblock are all the same.

        far_lower_diag[index] = -((n_mob(i,j) + n_mob(i+1, j))/2.)*Bn_negZ(i,j);
    }
}


//main lower diag
void Continuity_n::set_lower_diag()
{
    for (int index = 1; index <= num_elements-1; index++) {
        int i = index % N;
        int j = 1 + floor((index-1)/N);

        if (index % N == 0)
            lower_diag[index] = 0;      //probably don't need explicitely fill here, since auto intialized to 0
        else
            lower_diag[index] = -((n_mob(i,j) + n_mob(i,j+1))/2.)*Bn_negX(i,j);
    }
}


void Continuity_n::set_main_diag()
{
    for (int index = 1; index <= num_elements; index++) {
        int i = index % N;

        if (i == 0)
            i = N;

        int j = 1 + floor((index-1)/N);

        main_diag[index] = ((n_mob(i,j) + n_mob(i,j+1))/2.)*Bn_posX(i,j) + ((n_mob(i+1,j) + n_mob(i+1,j+1))/2.)*Bn_negX(i+1,j) + ((n_mob(i,j) + n_mob(i+1,j))/2.)*Bn_posZ(i,j) + ((n_mob(i,j+1) + n_mob(i+1,j+1))/2.)*Bn_negZ(i,j+1);
    }
}


void Continuity_n::set_upper_diag()
{
    for (int index = 1; index <= num_elements-1; index++) {
        int i = index % N;
        int j = 1 + floor((index-1)/N);

        if ((index-1 % N) == 0)
            upper_diag[index] = 0;
        else
            upper_diag[index] = -((n_mob(i+1,j) + n_mob(i+1,j+1))/2.)*Bn_posX(i+1,j);
    }
}


void Continuity_n::set_far_upper_diag()
{
    for (int index = 1; index <= num_elements-N; index++) {
        int i = index % N;

        if(i == 0)
            i = N;

        int j = 1 + floor((index-N)/N);

        far_upper_diag[index] = -((n_mob(i,j+1) + n_mob(i+1,j+1))/2.)*Bn_posZ(i,j+1);
    }
}











//---------------------------------------------------------------------------

void Continuity_n::set_rhs(const Eigen::MatrixXd &Un_matrix)
{
    for (int i = 1; i < rhs.size(); i++) {
        rhs[i] = -Cn*Un[i];
    }
    //BCs
    rhs[1] -= n_mob[0]*B_n2[1]*n_leftBC;
    rhs[rhs.size()-1] -= n_mob[rhs.size()]*B_n1[rhs.size()]*n_rightBC;
}

//------------------------
//Note: are using the V matrix for Bernoulli calculations.
//Makes it clearer to write indices in terms of (x,z) real coordinate values.

void Continuity_n::Bernoulli_n_X(const Eigen::MatrixXd &V_matrix)
{
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(num_cell+1,num_cell+1);

    for (int i = 1; i < num_cell+1; i++)             //Note: the indices are shifted by 1 from Matlab version (here bndry is at index 0)
        for (int j = 1; j < num_cell+1; j++)
            dV(i,j) =  V_matrix(i,j)-V_matrix(i-1,j);


    for (int i = 1; i < num_cell+1; i++) {           //note: the indexing done a bit different than Matlab (see 1D C++ version)
        for (int j = 1; j < num_cell+1; j++) {
            Bn_posX(i,j) = dV(i,j)/(exp(dV(i,j)) - 1.0);
            Bn_negX(i,j) = Bn_posX(i,j)*exp(dV(i,j));
        }
    }
}

void Continuity_n::Bernoulli_n_Z(const Eigen::MatrixXd &V_matrix)
{
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(num_cell+1,num_cell+1);

    for (int i = 1; i < num_cell+1; i++)
       for (int j = 1; j < num_cell+1; j++)
           dV(i,j) =  V_matrix(i,j)-V_matrix(i,j-1);

    for (int i = 1; i < num_cell+1; i++) {
        for (int j = 1; j < num_cell+1; j++) {

            Bn_posZ(i,j) = dV(i,j)/(exp(dV(i,j)) - 1.0);
            Bn_negZ(i,j) =  Bn_posZ(i,j)*exp(dV(i,j));
        }
    }
}
