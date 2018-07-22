#include "continuity_p.h"

Continuity_p::Continuity_p(const Parameters &params)
{
    num_elements = params.num_elements;
    num_cell_x = params.num_cell_x;
    num_cell_y = params.num_cell_y;
    num_cell_z = params.num_cell_z;
    Nx = num_cell_x - 1;
    Ny = num_cell_y - 1;
    Nz = num_cell_z - 1;

    p_matrix = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);  //note: this variable is updated every type find a new p_matrix inside main
    //later, can somehow make this more efficient

    //these diags vectors are not needed
//    main_diag.resize(num_elements+1);
//    upper_diag.resize(num_elements+1);
//    lower_diag.resize(num_elements+1);
//    far_lower_diag.resize(num_elements+1);
//    far_upper_diag.resize(num_elements+1);
    rhs.resize(num_elements+1);  //+1 b/c I am filling from index 1

    Bp_posX = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);
    Bp_negX = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);
    Bp_posY = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);
    Bp_negY = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);
    Bp_posZ = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);
    Bp_negZ = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);

    values  = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);
    values2 = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);
    values3 = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);
    values4 = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);

    Jp_Z = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);
    Jp_X = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);
    Jp_Y = Eigen::Tensor<double, 3> (num_cell_x+1, num_cell_y+1, num_cell_z+1);


    p_bottomBC.resize(num_cell_x+1, num_cell_y+1);
    p_topBC.resize(num_cell_x+1, num_cell_y+1);

    J_coeff_x = (q*Vt*params.N_dos*params.mobil)/params.dx;
    J_coeff_y = (q*Vt*params.N_dos*params.mobil)/params.dy;
    J_coeff_z = (q*Vt*params.N_dos*params.mobil)/params.dz;

    //------------------------------------------------------------------------------------------
    p_mob = Eigen::Tensor<double, 3> (num_cell_x+2, num_cell_y+2, num_cell_z+2);
    p_mob_avg_X = Eigen::Tensor<double, 3> (num_cell_x+2, num_cell_y+2, num_cell_z+2);
    p_mob_avg_Y = Eigen::Tensor<double, 3> (num_cell_x+2, num_cell_y+2, num_cell_z+2);
    p_mob_avg_Z = Eigen::Tensor<double, 3> (num_cell_x+2, num_cell_y+2, num_cell_z+2);
    p_mob.setConstant(params.p_mob_active/params.mobil);

    //Compute averaged mobilities
    for (int i = 1; i <= num_cell_x; i++) {
        for (int j = 1; j <= num_cell_y; j++) {
            for (int k = 1; k <= num_cell_z; k++) {
                p_mob_avg_X(i,j,k) = (p_mob(i,j,k) + p_mob(i,j+1,k) + p_mob(i,j,k+1) + p_mob(i,j+1,k+1))/4.;
                p_mob_avg_Y(i,j,k) = (p_mob(i,j,k) + p_mob(i+1,j,k) + p_mob(i,j,k+1) + p_mob(i+1,j,k+1))/4.;
                p_mob_avg_Z(i,j,k) = (p_mob(i,j,k) + p_mob(i+1,j,k) + p_mob(i,j+1,k) + p_mob(i+1,j+1,k))/4.;

                //add num_cell+2 values for i to account for extra bndry pt
                p_mob_avg_X(num_cell_x+1,j,k) = p_mob_avg_X(1,j,k);  //2 b/c the p_mob avg indices start from 2 (like Bernoullis)
                p_mob_avg_Y(num_cell_x+1,j,k) = p_mob_avg_Y(1,j,k);
                p_mob_avg_Z(num_cell_x+1,j,k) = p_mob_avg_Z(1,j,k);

                //add num_cell+2 values for j to account for extra bndry pt
                p_mob_avg_X(i,num_cell_y+1,k) = p_mob_avg_X(i,1,k);
                p_mob_avg_Y(i,num_cell_y+1,k) = p_mob_avg_Y(i,1,k);
                p_mob_avg_Z(i,num_cell_y+1,k) = p_mob_avg_Z(i,1,k);
            }
           //NOTE: use inside the device value for the num_cell+2 values b/c
           //these actually correspond to num_cell values, since with  Neuman
           //BC's, we use the value inside the device twice, when discretizing
           //at the boundary.
            p_mob_avg_X(i,j,num_cell_z+1) = p_mob(i,j,num_cell_z);        //assume the average p_mob at the bndry pt (top electrode) = same as p_mob just inside.
            p_mob_avg_Y(i,j,num_cell_z+1) = p_mob(i,j,num_cell_z);
            p_mob_avg_Z(i,j,num_cell_z+1) = p_mob(i,j,num_cell_z);
        }
    }

     //------------------------------------------------------------------------------------------
    Cp = (params.dx*params.dx)/(Vt*params.N_dos*params.mobil);  //can't use static, b/c dx wasn't defined as const, so at each initialization of Continuity_p object, new const will be made.

    //these BC's for now stay constant throughout simulation, so fill them once, upon Continuity_p object construction
    for (int j =  0; j <= num_cell_y; j++) {
        for (int i = 1; i <= num_cell_x+1; i++) {
            p_bottomBC(i,j) = 0;
            p_topBC(i,j) = 1;
        }
    }

    //allocate memory for the sparse matrix and rhs vector (Eig object)
    sp_matrix.resize(num_elements, num_elements);
    VecXd_rhs.resize(num_elements);   //only num_elements, b/c filling from index 0 (necessary for the sparse solver)

    //setup the triplet list for sparse matrix
     triplet_list.resize(11*num_elements);   //approximate the size that need
}

//------------------------------------------------------------------


//Calculates Bernoulli fnc values, then sets the diagonals and rhs
void Continuity_p::setup_eqn(const Eigen::Tensor<double, 3> &V_matrix, const std::vector<double> &Up, const Eigen::Tensor<double, 3> &p)
{
    trp_cnt = 0;  //reset triplet count
    Bernoulli_p(V_matrix);

    set_lowest_diag();
    set_lower_diag_Xs();   //lower diag corresponding to X direction finite differences
    set_lower_diag_Y_PBCs(); //lower diag corresponding to Y periodic boundary conditions
    set_lower_diag_Ys();
    set_main_lower_diag();
    set_main_diag();
    set_main_upper_diag();  //corresponds to Z direction finite differences
    set_upper_diag_Ys();
    set_upper_diag_Y_PBCs();
    set_upper_diag_Xs();
    set_highest_diag();

    set_rhs(Up);

    sp_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());    //sp_matrix is our sparse matrix
}

//------------------------------Setup Ap diagonals----------------------------------------------------------------
//X's left PBC
void Continuity_p::set_lowest_diag()
{
    values = -p_mob_avg_X*Bp_posX;

    int index = 1;
    int i = 1;     //since is PBC, this is always i = 1
    for (int j = 1; j <= Ny+1; j++) {
        for (int k = 1; k <= Nz; k++) {  //ONLY GOES TO Nz, b/c of Dirichlet BC's at top electrode (included in the matrix)..., all elements excep main diag need to be 0
            triplet_list[trp_cnt] = {index-1+(Nx)*(Nz+1)*(Ny+1), index-1, values(i,j,k)};  //note: don't need +1, b/c c++ values correspond directly to the inside pts
            //just  fill directly!! the triplet list. DON'T NEED THE DIAG VECTORS AT ALL!
            //RECALL, THAT the sparse matrices are indexed from 0 --> that's why have the -1's
            trp_cnt++;
            index = index +1;
        }
    }
}

//X's
void Continuity_p::set_lower_diag_Xs()
{
    values = -p_mob_avg_X*Bp_posX;

    int index = 1;
    for (int i = 2; i <= Nx+1; i++) {
        for (int j = 1; j <= Ny+1; j++) {
            for (int k = 1; k <= Nz; k++) {// only to Nz b/c of Dirichlet BCs
                triplet_list[trp_cnt] = {index-1+(Nz+1)*(Ny+1), index-1, values(i,j,k)};
                trp_cnt++;
                index = index +1;
            }
        }
    }
}


//Y's left PBCs
void Continuity_p::set_lower_diag_Y_PBCs()
{
    values = -p_mob_avg_Y*Bp_posY;

    int index = 1;
    int j = 1;   //always 1 b/c are bndry elements
    for (int i = 1; i <= Nx+1; i++) {
        for (int k = 1; k <= Nz; k++) { // only to Nz b/c of Dirichlet BCs
            triplet_list[trp_cnt] = {index-1+(Nz+1)*(Ny+1), index-1, values(i,j,k)};
            trp_cnt++;
            index = index +1;
        }
        index = index + (Nz+1)*(Ny);  // to skip the 0's subblocks
    }
}


//Y's
void Continuity_p::set_lower_diag_Ys()
{
    values = -p_mob_avg_Y*Bp_posY;

    int index = 1;
    for (int i = 1; i <= Nx+1; i++) {
        for (int j = 2; j <= Ny+1; j++) {
            for (int k = 1; k <= Nz; k++) {// only to Nz b/c of Dirichlet BCs
                triplet_list[trp_cnt] = {index-1+(Nz+1), index-1, values(i,j,k)};
                trp_cnt++;
                index = index +1;
            }
        }
        index = index + (Nz+1);  // to skip the 0's subblocks
    }
}



//main lower diag
void Continuity_p::set_main_lower_diag()
{
    values = -p_mob_avg_Z*Bp_posZ;

    int index = 1;
    for (int i = 1; i <= Nx+1; i++) {
        for (int j = 1; j <= Ny+1; j++) {
            for (int k = 2; k <= Nz; k++) {// only to Nz b/c of Dirichlet BCs
                triplet_list[trp_cnt] = {index, index-1, values(i,j,k)};
                trp_cnt++;
                index = index +1;
            }
            index = index + 1;  //skip the corner elements which are zero
        }
    }

}


//main diag
void Continuity_p::set_main_diag()
{
    values = p_mob_avg_Z*Bp_negZ + p_mob_avg_Y*Bp_negY + p_mob_avg_X*Bp_negX;
    values2 = p_mob_avg_X*Bp_posX;  //these have +1+1 in some elements, so need to be seperate
    values3 = p_mob_avg_Y*Bp_posY;
    values4 = p_mob_avg_Z*Bp_posZ;

    int index = 1;
    for (int i = 1; i <= Nx+1; i++) {
        for (int j = 1; j <= Ny+1; j++) {
            for (int k = 1; k <= Nz; k++) { // only to Nz b/c of Dirichlet BCs
                triplet_list[trp_cnt] = {index-1, index-1, values(i,j,k) + values2(i+1,j,k) + values3(i,j+1,k) + values4(i,j,k+1)};
                trp_cnt++;
                index = index +1;
            }
            //add the Dirichlet BC's element --> in matrix just have a 1
            int k = Nz+1;
            triplet_list[trp_cnt] = {index-1, index-1, 1};
            trp_cnt++;
            index = index + 1;
        }
    }
}


//main upper diag
void Continuity_p::set_main_upper_diag()
{
    values = -p_mob_avg_Z*Bp_negZ;

    int index = 1;  //note: unlike Matlab, can always start index at 1 here, b/c not using any spdiags fnc
    for (int i = 1; i <= Nx+1; i++) {
        for (int j = 1; j <= Ny+1; j++) {
            for (int k = 2; k <= Nz; k++) { // only to Nz b/c of Dirichlet BCs
                triplet_list[trp_cnt] = {index-1, index, values(i,j,k)};
                trp_cnt++;
                index = index +1;
            }
            index = index + 1; //to skip the 0 corner elements
        }
    }
}

//Y's
void Continuity_p::set_upper_diag_Ys()
{
    values = -p_mob_avg_Y*Bp_negY;

    int index = 1;
    for (int i = 1; i <= Nx+1; i++) {
        for (int j = 2; j <= Ny+1; j++) {
            for (int k = 1; k <= Nz; k++) { // only to Nz b/c of Dirichlet BCs
                triplet_list[trp_cnt] = {index-1, index-1+(Nz+1), values(i,j,k)};
                trp_cnt++;
                index = index +1;
            }
        }
        index = index + (Nz+1);  // to skip the 0's subblocks
    }
}


//Y right PBCs
void Continuity_p::set_upper_diag_Y_PBCs()
{
    values = -p_mob_avg_Y*Bp_negY;

    int index = 1;
    int j = Ny+1;  //corresponds to right y boundary
    for (int i = 1; i <= Nx+1; i++) {
        for (int k = 1; k <= Nz; k++) { // only to Nz b/c of Dirichlet BCs
            triplet_list[trp_cnt] = {index-1, index-1+(Nz+1), values(i,j,k)};
            trp_cnt++;
            index = index +1;
        }
        index = index + (Nz+1)*(Ny);  // to skip the 0's subblocks
    }
}


//X's
void Continuity_p::set_upper_diag_Xs()
{
    values = -p_mob_avg_X*Bp_negX;

    int index = 1;
    for (int i = 2; i <= Nx+1; i++) {
        for (int j = 1; j <= Ny+1; j++) {
            for (int k = 1; k <= Nz; k++) {// only to Nz b/c of Dirichlet BCs
                triplet_list[trp_cnt] = {index-1, index-1+(Nz+1)*(Ny+1), values(i,j,k)};
                trp_cnt++;
                index = index +1;
            }
        }
    }
}

//far upper diag X right PBC's
void Continuity_p::set_highest_diag()
{
    values = -p_mob_avg_X*Bp_negX;

    int index = 1;
    int i = Nx+1;     //corresponds to right boundary
    for (int j = 1; j <= Ny+1; j++) {
        for (int k = 1; k <= Nz; k++) {  // only to Nz b/c of Dirichlet BCs
            triplet_list[trp_cnt] = {index-1, index-1+(Nx)*(Nz+1)*(Ny+1), values(i,j,k)};
            trp_cnt++;
            index = index +1;
        }
    }
}

//---------------------------

void Continuity_p::Bernoulli_p(const Eigen::Tensor<double, 3> &V_matrix)
{
    Eigen::Tensor<double, 3> dV_X(num_cell_x+1, num_cell_y+1, num_cell_z+1),  dV_Y(num_cell_x+1, num_cell_y+1, num_cell_z+1),  dV_Z(num_cell_x+1, num_cell_y+1, num_cell_z+1);  //1 less than matlab, b/c my bndries are at 0 indices
    dV_X.setZero();
    dV_Y.setZero();
    dV_Z.setZero();

     for (int k = 1; k <= num_cell_z; k++) {
        for (int j = 1; j <= num_cell_y; j++) {
            for (int i = 1; i <= num_cell_x; i++) {
                dV_X(i,j,k) = V_matrix(i,j,k) - V_matrix(i-1,j,k);   //NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
                dV_Y(i,j,k) = V_matrix(i,j,k) - V_matrix(i,j-1,k);
                dV_Z(i,j,k) = V_matrix(i,j,k) - V_matrix(i,j,k-1);

                //add boundary case for num_cell+2 on Z value--> that dV = 0
                //from the Neuman boundary condition! And where we do have tip,
                //we assume tip voltage is same at +1 into tip, as at tip
                //contact bndry pt.
                dV_X(i,j, num_cell_z+1) = 0;
                dV_Y(i,j, num_cell_z+1) = 0;
                dV_Z(i,j, num_cell_z+1) = 0;

                //add on the wrap around dV's
                dV_X(i,num_cell_y+1,k) = V_matrix(i,0,k) - V_matrix(i-1,0,k);   //NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
                dV_Y(i,num_cell_y+1,k) = V_matrix(i,0,k) - V_matrix(i,num_cell_y,k);
                dV_Z(i,num_cell_y+1,k) = V_matrix(i,0,k) - V_matrix(i,0,k-1);

            }
            //add on the wrap around dV's and right boundary pt dV's
            dV_X(num_cell_x+1,j,k) = V_matrix(0,j,k) - V_matrix(num_cell_x,j,k);   //NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
            dV_Y(num_cell_x+1,j,k) = V_matrix(0,j,k) - V_matrix(0,j-1,k);  //use 0's b/c num_cell+2 is same as the 1st pt--> PBC's
            dV_Z(num_cell_x+1,j,k) = V_matrix(0,j,k) - V_matrix(0,j,k-1);
        }
        dV_X(num_cell_x+1,num_cell_y+1,k) = V_matrix(0,num_cell_y,k) - V_matrix(num_cell_x,num_cell_y,k);   //NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
        dV_Y(num_cell_x+1,num_cell_y+1,k) = V_matrix(0,num_cell_y,k) - V_matrix(0,num_cell_y-1,k);  //use 0's b/c num_cell+2 is same as the 1st pt--> PBC's
        dV_Z(num_cell_x+1,num_cell_y+1,k) = V_matrix(0,num_cell_y,k) - V_matrix(0,num_cell_y,k-1);
    }

    // matrix operations for Bernoulli's
    Bp_posX = dV_X/(dV_X.exp() - 1.0);  //pay attention!: is  . notation for exp() component wise operation
    Bp_negX = Bp_posX*dV_X.exp();

    Bp_posY = dV_Y/(dV_Y.exp() - 1.0);
    Bp_negY = Bp_posY*dV_Y.exp();

    Bp_posZ = dV_Z/(dV_Z.exp() - 1.0);
    Bp_negZ = Bp_posZ*dV_Y.exp();

    //I THINK THIS IS STILL OK TO USE--> only a ffefcts the small dV's--> i.e.
    //which are negligible. For larger dV's, we use the real Bernoulli eqn's...
    for (int i = 1; i <= num_cell_x+1; i++) {  //num_cell +1 for i and j b/c of the PBC's/including boundary pt--> but are setting to 1 anyway..., so doesn't matter much
        for (int j = 1; j <= num_cell_y+1; j++) {
            for (int k = 1; k <= num_cell_z+1; k++) {
                if(abs(dV_X(i,j,k)) < 1e-13) { //USING REAL TAYLOR EXPANSION HERE CAUSES BLOWUP!!!--> WHY?
                    Bp_posX(i,j,k) = 1;//1 -  dV_X(i,j,k)/2 + (dV_X(i,j,k))^2/12 - (dV_X(i,j,k))^4/720;
                    Bp_negX(i,j,k) =  1;//Bp_posX(i,j)*exp(dV_X(i,j,k));
                }
                if(abs(dV_Y(i,j,k)) < 1e-13) { //using real taylor expansion here works fine...
                    Bp_posY(i,j,k) = 1;//1 -  dV_Y(i,j,k)/2 + (dV_Y(i,j,k))^2/12 - (dV_Y(i,j,k))^4/720;
                    Bp_negY(i,j,k) = Bp_posY(i,j,k)*exp(dV_Y(i,j,k));
                }
                if(abs(dV_Z(i,j,k)) < 1e-13) { //using real taylor expansion here works fine...
                    Bp_posZ(i,j,k) = 1;//1 -  dV_Z(i,j,k)/2 + (dV_Z(i,j,k))^2/12 - (dV_Z(i,j,k))^4/720;
                    Bp_negZ(i,j,k) = Bp_posZ(i,j,k)*exp(dV_Z(i,j,k));
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
    for (int i = 1; i <= Nx+1; i++) {  //num_cell +1 for i and j b/c of the PBC's/including boundary pt--> but are setting to 1 anyway..., so doesn't matter much
        for (int j = 1; j <= Ny+1; j++) {
            for (int k = 1; k <= Nz+1; k++) {
                index++;
                if (k == 1) //bottom BC
                    rhs[index] += p_mob(i,j,0)*p_bottomBC(i,j)*Bp_posZ(i,j,k); //note: no +1's on Bp_posZ b/c I'm using shifted by 1 from Matlab, i.e. dV(1)-dV(0) = Bp(1)

                if (k == Nz + 1) //top BC
                    rhs[index] = p_topBC(i,j);

            }
        }
    }

    //set up VectorXd Eigen vector object for sparse solver
    for (int i = 1; i<=num_elements; i++) {
        VecXd_rhs(i-1) = rhs[i];   //fill VectorXd  rhs of the equation
    }
}


void Continuity_p::calculate_currents()
{
    for (int i = 1; i <= num_cell_x; i++) {
        for (int j = 1; j < num_cell_y; j++) {
            for (int k = 1; k < num_cell_z; k++) {
            Jp_Z(i,j,k) = -J_coeff_x * p_mob(i,j,k) * (p_matrix(i,j,k)*Bp_negZ(i,j,k) - p_matrix(i,j,k-1)*Bp_posZ(i,j,k));
            Jp_Y(i,j,k) = -J_coeff_y * p_mob(i,j,k) * (p_matrix(i,j,k)*Bp_negY(i,j,k) - p_matrix(i,j-1,k)*Bp_posY(i,j,k));
            Jp_X(i,j,k) = -J_coeff_z * p_mob(i,j,k) * (p_matrix(i,j,k)*Bp_negX(i,j,k) - p_matrix(i-1,j,k)*Bp_posX(i,j,k));
            }
        }
    }

}
