/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Solving 2D Poisson + Drift Diffusion semiconductor eqns for a solar cell using
%                      Scharfetter-Gummel discretization
%
%                         Written by Timofey Golubev
%
%     This includes the 2D poisson equation and 2D continuity/drift-diffusion
%     equations using Scharfetter-Gummel discretization. The Poisson equation
%     is solved first, and the solution of potential is used to calculate the
%     Bernoulli functions and solve the continuity eqn's.
%
%   Boundary conditions for Poisson equation are:
%
%     -a fixed voltage at (x,0) and (x, Nz) defined by V_bottomBC
%      and V_topBC which are defining the  electrodes
%
%    -insulating boundary conditions: V(0,z) = V(1,z) and
%     V(0,N+1) = V(1,N) (N is the last INTERIOR mesh point).
%     so the potential at the boundary is assumed to be the same as just inside
%     the boundary. Gradient of potential normal to these boundaries is 0.
%
%   Matrix equations are AV*V = bV, Ap*p = bp, and An*n = bn where AV, Ap, and An are sparse matrices
%   (generated using spdiag), for the Poisson and continuity equations.
%   V is the solution for electric potential, p is the solution for hole
%   density, n is solution for electron density
%   bV is the rhs of Poisson eqn which contains the charge densities and boundary conditions
%   bp is the rhs of hole continuity eqn which contains net generation rate
%   and BCs
%
%     The code as is will calculate data for a JV curve
%     as well as carrier densities, current densities, and electric field
%     distributions of a generic solar cell made of an active layer and electrodes.
%     More equations for carrier recombination can be easily added.f
%
%     Photogeneration rate will be inputed from gen_rate.inp file
%     (i.e. the output of an optical model can be used) or an analytic expression
%     for photogeneration rate can be added to photogeneration.cpp. Generation rate file
%     should contain num_cell-2 number of entries in a single column, corresponding to
%     the the generation rate at each mesh point (except the endpoints).
%
%     The code can also be applied to non-illuminated devices by
%     setting photogeneration rate to 0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>   //allows to use fill and min
#include <fstream>
#include <chrono>
#include <string>
#include <time.h>
#include <fstream>
#include <string>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include<Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include<Eigen/SparseLU>
#include <unsupported/Eigen/CXX11/Tensor>  //allows for 3D matrices (Tensors)

#include "constants.h"        //these contain physics constants only
#include "parameters.h"
#include "poisson.h"
#include "continuity_p.h"
//#include "continuity_n.h"
//#include "recombination.h"
//#include "photogeneration.h"
#include "Utilities.h"


int main()
{
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();  //start clock timer
    Parameters params;    //params is struct storing all parameters
    params.Initialize();  //reads parameters from file

    const int num_cell_x = params.num_cell_x;   //create a local num_cell so don't have to type params.num_cell everywhere
    const int num_cell_y = params.num_cell_y;
    const int num_cell_z = params.num_cell_z;

    const int num_V = static_cast<int>(floor((params.Va_max-params.Va_min)/params.increment))+1;  //floor returns double, explicitely cast to int
    params.tolerance_eq = params.tolerance_i;
    const int Nx = params.Nx;
    const int Ny = params.Ny;
    const int Nz = params.Nz;
    const int num_rows = (Nx+1)*(Ny+1)*(Nz+1);  //number of rows in the solution vectors (V, n, p)
    //NOTE: num_rows is the same as num_elements.
    //NOTE: we include the top BC inside the matrix and solution vectors to allow in future to use mixed BC's there.

    std::ofstream JV;
    JV.open("JV.txt");  //note: file will be created inside the build directory


    //-------------------------------------------------------------------------------------------------------
    //Initialize other vectors
    //Will use indicies for n and p... starting from 1 --> since is more natural--> corresponds to 1st node inside the device...
    //Note: these are Tensor type, so can use reshape on them, but these all are just a single column (i.e. the soln column from the matrix eqn).
    Eigen::Tensor<double, 3> n(num_rows+1, 1, 1), p(num_rows + 1, 1, 1), oldp(num_rows + 1, 1, 1), newp(num_rows + 1, 1, 1), oldn(num_rows + 1, 1, 1), newn(num_rows + 1, 1, 1);
    Eigen::Tensor<double, 3> oldV(num_rows + 1, 1, 1), newV(num_rows + 1, 1, 1), V(num_rows + 1, 1, 1);

    //create matrices to hold the V, n, and p values (including those at the boundaries) according to the (x,z) coordinates.
    //allows to write formulas in terms of coordinates
    Eigen::VectorXd soln_Xd(num_rows);  //vector for storing solutions to the  sparse solver (indexed from 0, so only num_rows size)

    //For the following, only need gen rate on insides, so N+1 size is enough
    std::vector<double> Un(num_rows+1); //will store generation rate as vector, for easy use in rhs
    std::vector<double> Up = Un;
    //Eigen::Tensor<double, 3> R_Langevin(N+1,N+1,N+1), PhotogenRate(N+1,N+1,N+1);
    Eigen::Tensor<double, 3> J_total_Z(num_cell_x+1, num_cell_y+1, num_cell_z+1), J_total_X(num_cell_x+1, num_cell_y+1, num_cell_z+1), J_total_Y(num_cell_x+1, num_cell_y+1, num_cell_z+1);                  //matrices for spacially dependent current
    Eigen::Tensor<double, 3> V_matrix(num_cell_x+1, num_cell_y+1, num_cell_z+1);  //Note: is actually a Tensor in Eigen
    Eigen::SparseMatrix<double> input; //for feeding input matrix into BiCGSTAB, b/c it crashes if try to call get matrix from the solve call.

    std::cout << Eigen::nbThreads( ) << std::endl;  //displays the # of threads that will be used by Eigen--> mine displays 8, but doesn't seem like it's using 8.

//------------------------------------------------------------------------------------
    //Construct objects
    Poisson poisson(params);
    //Recombo recombo(params);
    Continuity_p continuity_p(params);  //note this also sets up the constant top and bottom electrode BC's
    //Continuity_n continuity_n(params);  //note this also sets up the constant top and bottom electrode BC's
    //Photogeneration photogen(params, params.Photogen_scaling, params.GenRateFileName);
    Utilities utils;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::UpLoType::Lower, Eigen::AMDOrdering<int>> SCholesky; //Note using NaturalOrdering is much much slower

    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> SQR;
    Eigen::SparseLU<Eigen::SparseMatrix<double> >  poisson_LU, cont_n_LU, cont_p_LU;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> BiCGStab_solver;  //BiCGStab solver object

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::UpLoType::Lower|Eigen::UpLoType::Upper > cg;


//--------------------------------------------------------------------------------------------
    //Define boundary conditions and initial conditions. Note: electrodes are at the top and bottom.
    double Va = 0;
    poisson.set_V_bottomBC(params, Va);
    poisson.set_V_topBC(params, Va);

    //Initial conditions
    //std::vector<double> diff;
    //for (int x = 0; x <= num_cell; x++)
       //diff[x] = (poisson.get_V_topBC()[x] - poisson.get_V_bottomBC()[x])/num_cell;    //note, the difference can be different at different x values..., diff is in Z directiont

    //for now assume diff is constant everywhere...
    double diff = (poisson.get_V_topBC(0,0) - poisson.get_V_bottomBC(0,0))/num_cell_z;  //this is  calculated correctly

    for (int k = 1; k <= Nz+1; k++)
        for (int i = 1; i <= Nx+1; i++)
            for (int j = 1; j <= Ny+1; j++)
                V_matrix(i, j, k) = poisson.get_V_bottomBC(i+1,j+1) +  diff*k;

    //need to permute the matrix, to be consistent with the z,y,x ordering which I use for the matrices when solving.
    //NOTE: for Tensor shuffle to work, NEED TO EXPLICTELY CREATE AN ARRAY--> this isn't clear from the documentation
    Eigen::array<int, 3> permutation = {{2,1,0}}; //array should HAVE THE SPECIFIC type:  ptrdiff_t  ==> used for pointer arithmetic and array indexing.
    //ptrdiff_t is the signed integer type of the result of subtracting 2 pointers.

    V_matrix = V_matrix.shuffle(permutation);  //shuffle permuts the tensor. Note: dimensions are indexed from 0. This is supposed to swap x and z values...
    //Returns a copy of the input tensor whose dimensions have been reordered according to the specified permutation. The argument shuffle is an array of Index values. Its size is the rank of the input tensor. It must contain a permutation of 0, 1, ..., rank - 1.

    //reshape the matrix to a single column. //Note: even though reshaping to a column, V still must be a TENSOR type for this to work!
     Eigen::array<ptrdiff_t, 3> reshape_sizes = {{num_rows+1, 1, 1}};  //+1 b/c starts indexing from 0, but we use from 1
     V = V_matrix.reshape(reshape_sizes);  //Note: even though reshaping to a column, V still must be a TENSOR type for this to work!

    //Fill n and p with initial conditions (need for error calculation)
    double min_dense = 0;//continuity_n.get_n_bottomBC(1,1) < continuity_p.get_p_topBC(1,1) ? continuity_n.get_n_bottomBC(1,1):continuity_p.get_p_topBC(1,1);  //this should be same as std::min  fnc which doesn't work for some reason
    //double min_dense = std::min (continuity_n.get_n_bottomBC(1,1), continuity_p.get_p_topBC(1,1));  //Note: I defined the get fnc to take as arguments the i,j values... //Note: the bc's along bottom and top are currently uniform, so index, doesn't really matter.
    for (int i = 1; i<= num_rows; i++) {
        //n[i] = min_dense;
        p(i) = min_dense;  //NOTE: p is tensor now, so need () to access elements
    }

    //Convert the n and p to n_matrix and p_matrix
    //USE THE RESHAPE FUNCTION!
    Eigen::array<ptrdiff_t, 3> to_matrix_sizes{{Nx+1, Ny+1, Nz+1}};  //this is a reshape before solving, so keep in x,y,z format
    Eigen::Tensor<double, 3> p_matrix = p.reshape(to_matrix_sizes);
    continuity_p.set_p_matrix(p_matrix);  //save p_matrix to continuity_p member variable--> THIS MIGHT BE INEFFIIENT, BUT DO IT FOR NOW

    //////////////////////MAIN LOOP////////////////////////////////////////////////////////////////////////////////////////////////////////

    int iter, not_cnv_cnt, Va_cnt;
    bool not_converged;
    double error_np, old_error;  //this stores max value of the error and the value of max error from previous iteration
    std::vector<double> error_np_vector(num_rows+1);  //note: since n and p solutions are in vector form, can use vector form here also

    for (Va_cnt = 0; Va_cnt <= num_V +1; Va_cnt++) {  //+1 b/c 1st Va is the equil run
        not_converged = false;
        not_cnv_cnt = 0;
        if (params.tolerance > 1e-5)
            std::cerr<<"ERROR: Tolerance has been increased to > 1e-5" <<std::endl;

        if (Va_cnt==0) {
            params.use_tolerance_eq();  //relaxed tolerance for equil. run
            params.use_w_eq();
            Va = 0;
        }
        else {
            Va = params.Va_min+params.increment*(Va_cnt-1);
        }
        if (Va_cnt == 1) {
            params.use_tolerance_i();  //reset tolerance back
            params.use_w_i();
            //PhotogenRate = photogen.getPhotogenRate();    //otherwise PhotogenRate is pre-initialized to 0 in this main.cpp when declared
        }
        std::cout << "Va = " << Va <<std::endl;

        //Reset top and bottom BCs (outside of loop b/c don't change iter to iter)
        poisson.set_V_bottomBC(params, Va);
        poisson.set_V_topBC(params, Va);

        //-----------------------------------------------------------
        error_np = 1.0;
        iter = 0;

        while (error_np > params.tolerance) {
            //std::cout << "Va " << Va <<std::endl;

            //-----------------Solve Poisson Equation------------------------------------------------------------------
            poisson.set_rhs(p);  //this finds netcharge and sets rhs
            //std::cout << poisson.get_sp_matrix() << std::endl;
            oldV = V;

            if (iter == 0) { //INSTEAD OF HAVING IF here, can move these 2 lines, outside of the loop
                poisson_LU.analyzePattern(poisson.get_sp_matrix());  //by doing only on first iter, since pattern never changes, save a bit cpu
                poisson_LU.factorize(poisson.get_sp_matrix());
            }
            soln_Xd = poisson_LU.solve(poisson.get_rhs());


/*
            if (iter == 0) {  //This is slower than  LU
                SCholesky.analyzePattern(poisson.get_sp_matrix());
                SCholesky.factorize(poisson.get_sp_matrix());         //since numerical values of Poisson matrix don't change for 1 set of BC's, can factorize, just on 1st iter
            }
            soln_Xd = SCholesky.solve(poisson.get_rhs());
            */
            //std::cout << soln_Xd << std::endl;
             //std::cout << "Poisson solver error " << poisson.get_sp_matrix() * soln_Xd - poisson.get_rhs() << std::endl;

            //RECALL, I am starting my V vector from index of 1, corresponds to interior pts...
            for (int i = 1; i<=num_rows; i++) {
                newV(i,1,1) = soln_Xd(i-1);   //fill VectorXd  rhs of the equation
                //NOTE: newV is a Tensor now, so use () to access
            }

            //Mix old and new solutions for V
            if (iter > 0)
                V  = utils.linear_mix(params, newV, oldV);
            else
                V = newV;

            //reshape solution to a V_matrix
            V_matrix = V.reshape(to_matrix_sizes);

            //------------------------------Calculate Net Generation Rate----------------------------------------------------------

            //R_Langevin = recombo.ComputeR_Langevin(params,n,p);
            //FOR NOW CAN USE 0 FOR R_Langevin

            if (Va_cnt > 0) {
                for (int i = 1; i <= num_rows; i++) {
                    Up[i] = 0; //params.Photogen_scaling;  //This is what was used in Matlab version for testing.   photogen.getPhotogenRate()(i,j); //- R_Langevin(i,j);
                }
            }

            //--------------------------------Solve equation for p------------------------------------------------------------

            continuity_p.setup_eqn(poisson.get_V_matrix(), Up, p);
            //std::cout << continuity_p.get_sp_matrix() << std::endl;   //Note: get rhs, returns an Eigen VectorXd
            oldp = p;
/*
            input = continuity_p.get_sp_matrix();
            if (iter == 0)
                BiCGStab_solver.analyzePattern(input);
            BiCGStab_solver.factorize(input);  //this computes preconditioner, if use along with analyzePattern (for 1st iter)
            //BiCGStab_solver.compute(input);  //this computes the preconditioner..compute(input);
            soln_Xd = BiCGStab_solver.solve(continuity_p.get_rhs());
*/

            if (iter == 0 )
                cont_p_LU.analyzePattern(continuity_p.get_sp_matrix());
            cont_p_LU.factorize(continuity_p.get_sp_matrix());
            soln_Xd = cont_p_LU.solve(continuity_p.get_rhs());


            //save results back into n std::vector. RECALL, I am starting my V vector from index of 1, corresponds to interior pts...
            for (int i = 1; i<=num_rows; i++) {
                newp(i,1,1) = soln_Xd(i-1);   //newp is now a tensor....
            }

            //------------------------------------------------

            //if get negative p's or n's set them = 0
            for (int i = 1; i <= num_rows; i++) {
                if (newp(i,1,1) < 0.0) newp(i,1,1) = 0;
                //if (newn[i] < 0.0) newn[i] = 0;
            }

            //calculate the error
            old_error = error_np;

            //THIS CAN BE MOVED TO A FUNCTION IN UTILS
            for (int i = 1; i <= num_rows; i++) {
                if (newp(i,1,1)!=0) {
                    error_np_vector[i] = (abs(newp(i,1,1)-oldp(i,1,1)))/abs(oldp(i,1,1));
                }
            }
            error_np = *std::max_element(error_np_vector.begin()+1,error_np_vector.end());  //+1 b/c we are not using the 0th element
            std::fill(error_np_vector.begin(), error_np_vector.end(),0.0);  //refill with 0's so have fresh one for next iter

            //auto decrease w if not converging
            if (error_np >= old_error)
                not_cnv_cnt = not_cnv_cnt+1;
            if (not_cnv_cnt > 2000) {
                params.reduce_w();
                params.relax_tolerance();
                not_cnv_cnt = 0;
            }

            p = utils.linear_mix(params, newp, oldp);

            //note: top and bottom BC's don't need to be changed for now, since assumed to be constant... (they are set when initialize continuity objects)

            //Convert p to p_matrix
            p_matrix = p.reshape(to_matrix_sizes);
            continuity_p.set_p_matrix(p_matrix);  //update member variable

            //std::cout << error_np << std::endl;
            //std::cout << "weighting factor = " << params.w << std::endl << std::endl;

            iter = iter+1;
        }

        //-------------------Calculate Currents using Scharfetter-Gummel definition--------------------------

        //continuity_n.calculate_currents();
        continuity_p.calculate_currents();

        J_total_Z = continuity_p.get_Jp_Z();// + continuity_n.get_Jn_Z();
        J_total_X = continuity_p.get_Jp_X();// + continuity_n.get_Jn_X();
        J_total_Y = continuity_p.get_Jp_Y();// + continuity_n.get_Jn_Y();

        //---------------------Write to file----------------------------------------------------------------
        utils.write_details(params, Va, poisson.get_V_matrix(), p_matrix, J_total_Z, Up);  //note just write p twice for now, place holder for n
        if(Va_cnt >0) utils.write_JV(params, JV, iter, Va, J_total_Z);


    }//end of main loop

    JV.close();

    std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(finish-start);
    std::cout << "CPU time = " << time.count() << std::endl;

    return 0;
}
