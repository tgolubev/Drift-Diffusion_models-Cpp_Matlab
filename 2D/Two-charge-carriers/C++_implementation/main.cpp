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

#include <Eigen/Dense>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include<Eigen/SparseQR>
#include <Eigen/OrderingMethods>

#include "constants.h"        //these contain physics constants only
#include "parameters.h"
#include "poisson.h"
#include "continuity_p.h"
#include "continuity_n.h"
#include "recombination.h"
#include "photogeneration.h"
#include "thomas_tridiag_solve.h"
#include "Utilities.h"

int main()
{
    Parameters params;    //params is struct storing all parameters
    params.Initialize();  //reads parameters from file

    const int num_cell = params.num_cell;   //create a local num_cell so don't have to type params.num_cell everywhere

    const int num_V = static_cast<int>(floor((params.Va_max-params.Va_min)/params.increment))+1;  //floor returns double, explicitely cast to int
    params.tolerance_eq = 100.*params.tolerance_i;
    const int N = params.num_cell -1;
    const int num_rows = N*N;  //number of rows in the solution vectors (V, n, p)

    std::ofstream JV;
    JV.open("JV.txt");  //note: file will be created inside the build directory

    //-------------------------------------------------------------------------------------------------------
    //Initialize other vectors
    //Will use indicies for n and p... starting from 1 --> since is more natural--> corresponds to 1st node inside the device...
    //USE SIZE num_rows + 1 b/c don't want to use the 0th index. Mathematically convenient to index the soln vectors from 1 to num_rows
    //NOTE: ALL THESE INCLUDE THE INTERIOR ELEMENTS ONLY
    std::vector<double> n(num_rows+ 1), p(num_rows+ 1), oldp(num_rows+ 1), newp(num_rows), oldn(num_rows+ 1), newn(num_rows+ 1);
    std::vector<double> oldV(num_rows+ 1), newV(num_rows+ 1), V(num_rows+ 1);

    //std::vector<double> Jp(num_rows),Jn(num_cell), J_total(num_cell);

    //create matrices to hold the V, n, and p values (including those at the boundaries) according to the (x,z) coordinates.
    //Allows to write formulas in terms of coordinates.
    //Note: these matrices are indexed from 0 (corresponds to the BC).
    //last index I need is num_cell - 1 = N
    Eigen::MatrixXd V_matrix = Eigen::MatrixXd::Zero(N+1,N+1);
    //Eigen::MatrixXd n_matrix = Eigen::MatrixXd::Zero(N+1,N+1);
    //Eigen::MatrixXd p_matrix = Eigen::MatrixXd::Zero(N+1,N+1);
    Eigen::MatrixXd Un_matrix = Eigen::MatrixXd::Zero(N+1,N+1);
    Eigen::MatrixXd Up_matrix = Eigen::MatrixXd::Zero(N+1,N+1);
    Eigen::MatrixXd R_Langevin(N+1,N+1), PhotogenRate(N+1,N+1);  //store the results of these..

    //define initial conditions as min value of BCs
    //double min_dense = std::min(continuity_n.get_n_bottomBC()[1], continuity_p.get_p_topBC()[1]);  //just use 1st index, the BC's are uniform for IC

    //initially we make the n and p densities inside the device be the same, so netcharge = 0
    Eigen::MatrixXd netcharge = Eigen::MatrixXd::Zero(num_cell+1,num_cell+1);

//------------------------------------------------------------------------------------
    //Construct objects
    Poisson poisson(params, netcharge);
    Recombo recombo(params);
    Continuity_p continuity_p(params);  //note this also sets up the constant top and bottom electrode BC's
    Continuity_n continuity_n(params);  //note this also sets up the constant top and bottom electrode BC's
    Photogeneration photogen(params, params.Photogen_scaling, params.GenRateFileName);
    Utilities utils;

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
    double diff = (poisson.get_V_topBC()[0] - poisson.get_V_bottomBC()[0])/num_cell;
    int index = 0;
    for (int j = 1; j <= N; j++) {//  %corresponds to z coord
        index++;
        V[index] = diff*j;
        for (int i = 2; i <= N; i++) {//  %elements along the x direction assumed to have same V
            index++;
            V[index] = V[index-1];
        }
    }

    //side BCs, insulating BC's
    poisson.set_V_leftBC(V);
    poisson.set_V_rightBC(V);

    //-----------------------

    //NOTE: the n and p top and bottom BC's are set when initialize a continuity_n and continuity_p object

    poisson.setup_matrix();  //outside of loop since matrix never changes

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();  //start clock timer

    //////////////////////MAIN LOOP////////////////////////////////////////////////////////////////////////////////////////////////////////

    int iter, not_cnv_cnt, Va_cnt;
    bool not_converged;
    double error_np, old_error;  //this stores max value of the error and the value of max error from previous iteration
    std::vector<double> error_np_vector(num_cell+1);  //note: since n and p solutions are in vector form, can use vector form here also

    for (Va_cnt = 0; Va_cnt <= num_V +1; Va_cnt++) {  //+1 b/c 1st Va is the equil run
        not_converged = false;
        not_cnv_cnt = 0;
        if (params.tolerance > 1e-5) {
            std::cerr<<"ERROR: Tolerance has been increased to > 1e-5" <<std::endl;
        }
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
            PhotogenRate = photogen.getPhotogenRate();    //otherwise PhotogenRate is pre-initialized to 0 in this main.cpp when declared
        }
        std::cout << "Va = " << Va <<std::endl;

        //Apply the voltage boundary conditions
        poisson.set_V_bottomBC(params, Va);
        poisson.set_V_topBC(params, Va);

        error_np = 1.0;
        iter = 0;
        while (error_np > params.tolerance) {
            std::cout << "error np " << error_np <<std::endl;
            std::cout << "Va " << Va <<std::endl;

            //-----------------Solve Poisson Equation------------------------------------------------------------------

            poisson.set_rhs(netcharge);
            oldV = V;

            //MATRIX SOVE WITH EIGEN

            //Mix old and new solutions for V
            if (iter > 0)
                V  = utils.linear_mix(params, newV, oldV);
            else
                V = newV;


            //conversion from V vector to V_matrix (holds only the inside values of V: indices 1 through N
            for (int index = 1; index <= num_rows; index++) {
                int i = index % N;    //this gives the i value for matrix
                if (i == 0) i = N;

                int j = 1 + static_cast<int>(floor((index-1)/N));  // j value for matrix

                V_matrix(i,j) = V[index];
            }

            //------------------------------Calculate Net Generation Rate----------------------------------------------------------

            //R_Langevin = recombo.ComputeR_Langevin(params,n,p);

            //------------//////////////////////////////////////////////////////////////
            //NEED AN R_Langevin in Matrix form--> do the conversion within the function....
            //////////////////////////////////////////////////////////////////


            //FOR NOW CAN USE 0 FOR R_Langevin

            for (int i = 1; i <= N; i++)
                for (int j = 1; j <= N; i++)
                    Un_matrix(i,j) = photogen.getPhotogenRate()(i,j); //- R_Langevin(i,j);

            Up_matrix = Un_matrix;


            //--------------------------------Solve equations for n and p------------------------------------------------------------ 

            continuity_n.setup_eqn(V_matrix, Un_matrix);
            oldn = n;
            //MATRIX SOVE WITH EIGEN

            continuity_p.setup_eqn(V_matrix, Up_matrix);
            oldp = p;
           //MATRIX SOVE WITH EIGEN

            //if get negative p's or n's set them = 0
            for (int i = 1; i < num_cell; i++) {
                if (newp[i] < 0.0) newp[i] = 0;
                if (newn[i] < 0.0) newn[i] = 0;
            }

            //calculate the error
            old_error = error_np;
            for (int i = 1; i < num_cell; i++) {
                if (newp[i]!=0 && newn[i] !=0) {
                    error_np_vector[i] = (abs(newp[i]-oldp[i]) + abs(newn[i]-oldn[i]))/abs(oldp[i]+oldn[i]);
                }
            }
            error_np = *std::max_element(error_np_vector.begin(),error_np_vector.end());
            std::fill(error_np_vector.begin(), error_np_vector.end(),0.0);  //refill with 0's so have fresh one for next iter

            //auto decrease w if not converging
            if (error_np >= old_error)
                not_cnv_cnt = not_cnv_cnt+1;
            if (not_cnv_cnt > 2000) {
                //params.reduce_w();
                params.relax_tolerance();
                not_cnv_cnt = 0;
            }

            p = utils.linear_mix(params, newp, oldp);
            n = utils.linear_mix(params, newn, oldn);

            //Apply continuity equation  BC's
            //set the  side BC's
            continuity_n.set_n_leftBC(n);
            continuity_n.set_n_rightBC(n);
            continuity_p.set_p_leftBC(p);
            continuity_p.set_p_rightBC(p);
            //note: top and bottom BC's don't need to be changed for now, since assumed to be constant... (they are set when initialize continuity objects)

            iter = iter+1;
        }

        std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(finish-start);
        std::cout << "1 Va CPU time = " << time.count() << std::endl;

        //-------------------Calculate Currents using Scharfetter-Gummel definition--------------------------
        /*
        for (int i = 1; i < num_cell; i++) {
            Jp[i] = -(q*Vt*params.N*params.mobil/params.dx) * continuity_p.get_p_mob()[i] * (p[i]*continuity_p.get_B_p2()[i] - p[i-1]*continuity_p.get_B_p1()[i]);
            Jn[i] =  (q*Vt*params.N*params.mobil/params.dx) * continuity_n.get_n_mob()[i] * (n[i]*continuity_n.get_B_n1()[i] - n[i-1]*continuity_n.get_B_n2()[i]);
            J_total[i] = Jp[i] + Jn[i];
        }
        */

        //---------------------Write to file----------------------------------------------------------------
        //utils.write_details(params, Va, V, p, n, J_total, Un, PhotogenRate, R_Langevin);
        //if(Va_cnt >0) utils.write_JV(params, JV, iter, Va, J_total);

    }//end of main loop
    JV.close();

    return 0;
}
