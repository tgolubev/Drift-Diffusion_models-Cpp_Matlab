#include "continuity_n.h"

Continuity_n::Continuity_n(const Parameters &params)
{
    num_elements = params.num_elements; //note: num_elements is same thing as num_rows in main.cpp
    N = params.num_cell - 1;
    num_cell = params.num_cell;

    main_diag.resize(num_elements+1);
    upper_diag.resize(num_elements);
    lower_diag.resize(num_elements);
    far_lower_diag.resize(num_elements-N+1);
    far_upper_diag.resize(num_elements-N+1);
    rhs.resize(num_elements+1);  //+1 b/c I am filling from index 1

   n_bottomBC.resize(num_cell+1);
   n_topBC.resize(num_cell+1);
   n_leftBC.resize(num_cell+1);
   n_rightBC.resize(num_cell+1);

   Bn_posX.resize(num_cell+1, num_cell+1);  //allocate memory for the matrix object
   Bn_negX.resize(num_cell+1, num_cell+1);
   Bn_posZ.resize(num_cell+1, num_cell+1);  //allocate memory for the matrix object
   Bn_negZ.resize(num_cell+1, num_cell+1);

   // //MUST FILL WITH THE VALUES OF n_mob!!  WILL NEED TO MODIFY THIS WHEN HAVE SPACE VARYING
   n_mob = (params.n_mob_active/params.mobil)*Eigen::MatrixXd::Ones(num_cell+1, num_cell+1);

   Cn = params.dx*params.dx/(Vt*params.N_dos*params.mobil);

   //these BC's for now stay constant throughout simulation, so fill them once, upon Continuity_n object construction
   for (int j =  0; j <= num_cell; j++) {
      n_bottomBC[j] = params.N_LUMO*exp(-(params.E_gap-params.phi_a)/Vt)/params.N_dos;
      n_topBC[j] = params.N_LUMO*exp(-params.phi_c/Vt)/params.N_dos;
   }

   //allocate memory for the sparse matrix and rhs vector (Eig object)
   sp_matrix.resize(num_elements, num_elements);
   VecXd_rhs.resize(num_elements);   //only num_elements, b/c filling from index 0 (necessary for the sparse solver)

   //setup the triplet list for sparse matrix
    //triplet_list.resize(5*num_elements);   //approximate the size that need         // list of non-zeros coefficients in triplet form(row index, column index, value)


}

//----------------------------------------------------------
//Set BC's

void Continuity_n::set_n_leftBC(const std::vector<double> &n)
{
    for (int j = 1; j <= N; j++) {
        n_leftBC[j] = n[(j-1)*N + 1];
    }


}

void Continuity_n::set_n_rightBC(const std::vector<double> &n)
{
    for (int j = 1; j <= N; j++) {
         n_rightBC[j]= n[j*N];
    }
}


//Calculates Bernoulli fnc values, then sets the diagonals and rhs
//use the V_matrix for setup, to be able to write equations in terms of (x,z) coordingates
void Continuity_n::setup_eqn(const Eigen::MatrixXd &V_matrix, const Eigen::MatrixXd &Un_matrix, const std::vector<double> &n)
{
    //trp_cnt = 0;  //reset triplet count --> used for setting up Eigen triplets
    Bernoulli_n_X(V_matrix);
    Bernoulli_n_Z(V_matrix);
    set_far_lower_diag();
    set_lower_diag();
    set_main_diag();
    set_upper_diag();
    set_far_upper_diag();
    set_n_rightBC(n);
    set_n_leftBC(n);
    set_rhs(Un_matrix);


    /////////////////////////////
    //NOTE: THIS CAN BE MOVED INTO THE SET DIAGONALS FUNCTIONS.... could be bit more efficient and clean
    /////////////////////////////
    //generate triplets for Eigen sparse matrix
    //setup the triplet list for sparse matrix
    typedef Eigen::Triplet<double> Trp;
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

//-------------------------------Setup An diagonals (Continuity/drift-diffusion solve)-----------------------------

void Continuity_n::set_far_lower_diag()
{
    //Lowest diagonal: corresponds to V(i, j-1)
    for (int index = 1; index <=N*(N-1); index++) {      //(1st element corresponds to Nth row  (number of elements = N*(N-1)
        int i = index % N;

        if (i ==0)                //the multiples of N correspond to last index
            i = N;

        int j = 1 + static_cast<int>(floor((index-1)/N));    //this is the y index of V which element corresponds to. 1+ floor(index/4)determines which subblock this corresponds to and thus determines j, since the j's for each subblock are all the same.

        far_lower_diag[index] = -((n_mob(i,j) + n_mob(i+1, j))/2.)*Bn_negZ(i,j);

        //triplet_list[trp_cnt] = {i-1, j-1, far_lower_diag[index]};  //i-1 and j-1 b/c Eigen matrices must be indexed from 0
        //trp_cnt++;
    }
}


//main lower diag
void Continuity_n::set_lower_diag()
{
    for (int index = 1; index <= num_elements-1; index++) {
        int i = index % N;
        int j = 1 + static_cast<int>(floor((index-1)/N));

        if (index % N == 0)
            lower_diag[index] = 0;      //probably don't need explicitely fill here, since auto intialized to 0
        else
            lower_diag[index] = -((n_mob(i,j) + n_mob(i,j+1))/2.)*Bn_negX(i,j);

        //triplet_list[trp_cnt] = {i-1, j-1, lower_diag[index]};
        //trp_cnt++;
    }
}


void Continuity_n::set_main_diag()
{
    for (int index = 1; index <= num_elements; index++) {
        int i = index % N;

        if (i == 0)
            i = N;

        int j = 1 + static_cast<int>(floor((index-1)/N));

        main_diag[index] = ((n_mob(i,j) + n_mob(i,j+1))/2.)*Bn_posX(i,j) + ((n_mob(i+1,j) + n_mob(i+1,j+1))/2.)*Bn_negX(i+1,j) + ((n_mob(i,j) + n_mob(i+1,j))/2.)*Bn_posZ(i,j) + ((n_mob(i,j+1) + n_mob(i+1,j+1))/2.)*Bn_negZ(i,j+1);

        //triplet_list[trp_cnt] = {i-1, j-1, main_diag[index]};
        //trp_cnt++;
    }
}


void Continuity_n::set_upper_diag()
{
    for (int index = 1; index <= num_elements-1; index++) {
        int i = index % N;
        int j = 1 + static_cast<int>(floor((index-1)/N));

        if ((index % N) == 0)
            upper_diag[index] = 0;
        else
            upper_diag[index] = -((n_mob(i+1,j) + n_mob(i+1,j+1))/2.)*Bn_posX(i+1,j);

        //triplet_list[trp_cnt] = {i-1, j-1, upper_diag[index]};
        //trp_cnt++;
    }


}


void Continuity_n::set_far_upper_diag()
{
    for (int index = 1; index <= num_elements-N; index++) {
        int i = index % N;

        if(i == 0)
            i = N;

        int j = 1 + static_cast<int>(floor((index-N)/N));

        far_upper_diag[index] = -((n_mob(i,j+1) + n_mob(i+1,j+1))/2.)*Bn_posZ(i,j+1);

        //triplet_list[trp_cnt] = {i-1, j-1, far_upper_diag[index]};
        //trp_cnt++;
    }
}

//---------------------------------------------------------------------------

void Continuity_n::set_rhs(const Eigen::MatrixXd &Un_matrix)
{
    int index = 0;

    for (int j = 1; j <= N; j++) {
        if (j ==1)  {//different for 1st subblock
            for (int i = 1; i <= N; i++) {
                index++;
                if (i==1) {     //1st element has 2 BC's
                    rhs[index] = Cn*Un_matrix(i,j) + n_mob(i,j)*(Bn_negX(i,j)*n_leftBC[1] + Bn_negZ(i,j)*n_bottomBC[i]);  //NOTE: rhs is +Cp*Un_matrix, b/c diagonal elements are + here, flipped sign from 1D version

                }
                else if (i==N)
                    rhs[index] = Cn*Un_matrix(i,j) + n_mob(i,j)*(Bn_negZ(i,j)*n_bottomBC[i] + Bn_posX(i+1,j)*n_rightBC[1]);
                else
                    rhs[index] = Cn*Un_matrix(i,j) + n_mob(i,j)*Bn_negZ(i,j)*n_bottomBC[i];
            }
        } else if (j == N) {      //different for last subblock
            for (int i = 1; i <= N; i++) {
                index++;
                if (i==1)  //1st element has 2 BC's
                    rhs[index] = Cn*Un_matrix(i,j) + n_mob(i,j)*(Bn_negX(i,j)*n_leftBC[N] + Bn_posZ(i,j+1)*n_topBC[i]);
                else if (i==N)
                        rhs[index] = Cn*Un_matrix(i,j) + n_mob(i,j)*(Bn_posX(i+1,j)*n_rightBC[N] + Bn_posZ(i,j+1)*n_topBC[i]);
                else
                rhs[index] = Cn*Un_matrix(i,j) + n_mob(i,j)*Bn_posZ(i,j+1)*n_topBC[i];
            }
        } else {     //interior subblocks
            for (int i = 1; i <= N; i++) {
                index++;
                if(i==1)
                    rhs[index] = Cn*Un_matrix(i,j) + n_mob(i,j)*Bn_negX(i,j)*n_leftBC[j];
                else if(i==N)
                        rhs[index] = Cn*Un_matrix(i,j) + n_mob(i,j)*Bn_posX(i+1,j)*n_rightBC[j];
                else
                rhs[index] = Cn*Un_matrix(i,j);
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

void Continuity_n::Bernoulli_n_X(const Eigen::MatrixXd &V_matrix)
{
    Eigen::MatrixXd dV = Eigen::MatrixXd::Zero(num_cell+1,num_cell+1);

    for (int i = 1; i < num_cell+1; i++)             //Note: the indices are shifted by 1 from Matlab version (here bndry is at index 0)
        for (int j = 1; j < num_cell+1; j++)
            dV(i,j) =  V_matrix(i,j)-V_matrix(i-1,j);


    for (int i = 1; i < num_cell+1; i++) {           //note: the indexing done a bit different than Matlab (see 1D C++ version)
        for (int j = 1; j < num_cell+1; j++) {
             if (abs(dV(i,j)) < 1e-13) {        //to prevent blowup due  to 0 denominator
                 Bn_posX(i,j) = 1 - dV(i,j)/2. + (dV(i,j)*dV(i,j))/12. - pow(dV(i,j), 4)/720.;
                 Bn_negX(i,j) =  Bn_posX(i,j)*exp(dV(i));
             } else {
                Bn_posX(i,j) = dV(i,j)/(exp(dV(i,j)) - 1.0);
                Bn_negX(i,j) = Bn_posX(i,j)*exp(dV(i,j));
             }
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
            if (abs(dV(i,j)) < 1e-13) {        //to prevent blowup due  to 0 denominator
                Bn_posZ(i,j) = 1 - dV(i,j)/2. + (dV(i,j)*dV(i,j))/12. - pow(dV(i,j), 4)/720.;
                Bn_negZ(i,j) =  Bn_posZ(i,j)*exp(dV(i));
            } else {
               Bn_posZ(i,j) = dV(i,j)/(exp(dV(i,j)) - 1.0);
               Bn_negZ(i,j) = Bn_posZ(i,j)*exp(dV(i,j));
            }
        }
    }
}
