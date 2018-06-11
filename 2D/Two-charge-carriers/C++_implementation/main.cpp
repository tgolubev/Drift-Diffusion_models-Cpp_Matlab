/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Solving 1D Poisson + Drift Diffusion semiconductor eqns for a solar cell using
%                      Scharfetter-Gummel discretization
%
%                         Written by Timofey Golubev
%                     Version 2.0  6/7/2018 (object oriented + file input)
%                     Version 1.0   5/26/17  (no objects)
%
%
%     The code as is will calculate and plot a JV curve
%     as well as carrier densities, current densities, and electric field
%     distributions of a generic solar cell made of an active layer and electrodes.
%     More equations for carrier recombination can be easily added.
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

#include "constants.h"        //these contain physics constants only
#include "parameters.h"
#include "poisson.h"
#include "continuity_p.h"
#include "continuity_n.h"
#include "recombination.h"
#include "photogeneration.h"
#include "thomas_tridiag_solve.h"
#include "bernoulli.h"


int main()
{
    Parameters params;    //params is struct storing all parameters
    params.Initialize();  //reads parameters from file

    const int num_cell = params.num_cell;   //create a local num_cell so don't have to type params.num_cell everywhere
    double Vbi = params.WF_anode - params.WF_cathode +params.phi_a +params.phi_c;
    int num_V = static_cast<int>(floor((params.Va_max-params.Va_min)/params.increment))+1;  //floor returns double, explicitely cast to int
    params.tolerance_eq = 100.*params.tolerance_i;

    std::ofstream JV, VaData;
    JV.open("JV.txt");  //note: file will be created inside the build directory

    //-------------------------------------------------------------------------------------------------------
//CAN MOVE THIS TO INITIALIZE

    //Fill the RELATIVE dielectric constant vector (can be position dependent, as long as piecewise constant)
    std::vector<double> epsilon(num_cell+1);
    std::fill(epsilon.begin(), epsilon.end(), params.eps_active);

    //Fill the mobilities vectors (can be position dependent, as long as piecewise constant)
    std::vector<double> p_mob(num_cell+1), n_mob(num_cell+1);
    std::fill(p_mob.begin(), p_mob.end(), params.p_mob_active);
    std::fill(n_mob.begin(), n_mob.end(), params.n_mob_active);
    for(int i = 0;i<= num_cell;i++){
        p_mob[i] = p_mob[i]/params.mobil;
        n_mob[i] = n_mob[i]/params.mobil;
    }

    //-------------------------------------------------------------------------------------------------------

    //Initialize other vectors
    //Will use indicies for n and p... starting from 1 --> since is more natural--> corresponds to 1st node inside the device...
    std::vector<double> n(num_cell), p(num_cell), oldp(num_cell), newp(num_cell), oldn(num_cell), newn(num_cell), oldV(num_cell+1), newV(num_cell+1), V(num_cell+1);
    std::vector<double> Un(num_cell), Up(num_cell), R_Langevin(num_cell), PhotogenRate(num_cell);  //store the results of these..
    std::vector<double> Jp(num_cell),Jn(num_cell), J_total(num_cell);
    std::vector<double> error_np_vector(num_cell);  //for storing errors between iterations
    std::vector<double> B_n1(num_cell+1), B_n2(num_cell+1), B_p1(num_cell+1), B_p2(num_cell+1); //vectors for storing Bernoulli fnc values

    //Construct objects
    Poisson poisson(params);
    Recombo recombo(params, params.k_rec_input);
    Continuity_p continuity_p(params);
    Continuity_n continuity_n(params);
    Photogeneration photogen(params, params.Photogen_scaling, params.GenRateFileName);

    //Initial conditions
    double min_dense = std::min(continuity_n.get_n_leftBC(),  continuity_p.get_p_rightBC());
    std::fill(n.begin()+1, n.end(), min_dense);
    std::fill(p.begin()+1, p.end(), min_dense);

    double V_leftBC = -((Vbi)/(2*Vt)-params.phi_a/Vt);
    double V_rightBC = (Vbi)/(2*Vt)-params.phi_c/Vt;
    double diff = (V_rightBC-V_leftBC)/num_cell;
    V[0] = V_leftBC;  //fill V(0) here for use in Beroulli later
    for(int i = 1; i<=num_cell-1; i++){
        V[i] = V[i-1] + diff;
    }
    V[num_cell] = V_rightBC;

    poisson.setup_matrix(epsilon);  //outside of loop since matrix never changes

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();  //start clock timer

    //////////////////////MAIN LOOP////////////////////////////////////////////////////////////////////////////////////////////////////////

    int iter, not_cnv_cnt, Va_cnt;
    double error_np, old_error;
    bool not_converged;
    double Va;

    for(Va_cnt = 0; Va_cnt <=num_V +1;Va_cnt++){  //+1 b/c 1st Va is the equil run
        not_converged = false;
        not_cnv_cnt = 0;
        if(params.tolerance > 1e-5){
            std::cerr<<"ERROR: Tolerance has been increased to > 1e-5" <<std::endl;
        }
        if(Va_cnt==0){
            params.use_tolerance_eq();  //relaxed tolerance for equil. run
            params.use_w_eq();
            Va = 0;
        }
        else{
            Va = params.Va_min+params.increment*(Va_cnt-1);
        }
        if(Va_cnt ==1){
            params.use_tolerance_i();  //reset tolerance back
            params.use_w_i();
            PhotogenRate = photogen.getPhotogenRate();    //otherwise PhotogenRate is pre-initialized to 0 in this main.cpp when declared
        }
        std::cout << "Va = " << Va <<std::endl;

        //Apply the voltage boundary conditions
        V_leftBC = -((Vbi  -Va)/(2*Vt)-params.phi_a/Vt);
        V_rightBC = (Vbi- Va)/(2*Vt) - params.phi_c/Vt;
        V[0] = V_leftBC;
        V[num_cell] = V_rightBC;

        error_np = 1.0;
        iter = 0;
        while(error_np > params.tolerance){
            std::cout << "error np " << error_np <<std::endl;
            std::cout << "Va " << Va <<std::endl;

            //-----------------Solve Poisson Equation------------------------------------------------------------------

            poisson.set_rhs(epsilon, n, p, V_leftBC, V_rightBC);
            oldV = V;
            newV = Thomas_solve(poisson.get_main_diag(), poisson.get_upper_diag(), poisson.get_lower_diag(), poisson.get_rhs());
            //add on the BC's --> b/c matrix solver just outputs the insides...
            newV[0] = V[0];
            newV[num_cell] = V[num_cell];

            //Mix old and new solutions for V
            if(iter>0){
                for(int i =  1;i<num_cell;i++){
                    V[i] = newV[i]*params.w + oldV[i]*(1.-params.w);
                }
            }
            else V = newV;

            //------------------------------Calculate Net Generation Rate----------------------------------------------------------

            R_Langevin = recombo.ComputeR_Langevin(params,n,p);

            for(int i = 1;i<num_cell;i++){
                Un[i] = PhotogenRate[i] - R_Langevin[i];
            }
            Up = Un;

            //--------------------------------Solve equations for n and p------------------------------------------------------------ 

            //setup the matrices
            BernoulliFnc_n(V, B_n1, B_n2);  //fills B_p1 and B_p2 with the values of Bernoulli fnc
            continuity_n.setup_eqn(n_mob, B_n1, B_n2, Un);
            oldn = n;
            newn = Thomas_solve(continuity_n.get_main_diag(), continuity_n.get_upper_diag(), continuity_n.get_lower_diag(), continuity_n.get_rhs());

            BernoulliFnc_p(V, B_p1, B_p2);
            continuity_p.setup_eqn(p_mob, B_p1, B_p2, Up);
            oldp = p;
            newp = Thomas_solve(continuity_p.get_main_diag(), continuity_p.get_upper_diag(), continuity_p.get_lower_diag(), continuity_p.get_rhs());

            //if get negative p's or n's set them = 0
            for(int i = 1; i<num_cell;i++){
                if(newp[i]<0.0) newp[i] = 0;
                if(newn[i]<0.0) newn[i] = 0;
            }

            //calculate the error
            old_error = error_np;
            for (int i = 1;i<num_cell;i++){
                if(newp[i]!= 0 && newn[i] !=0){
                    error_np_vector[i] = (abs(newp[i]-oldp[i])+abs(newn[i]-oldn[i]))/abs(oldp[i]+oldn[i]);
                }
            }
            error_np = *std::max_element(error_np_vector.begin(),error_np_vector.end());
            std::fill(error_np_vector.begin(), error_np_vector.end(),0.0);  //refill with 0's so have fresh one for next iter

            //auto decrease w if not converging
            if(error_np>=old_error) not_cnv_cnt = not_cnv_cnt+1;
            if(not_cnv_cnt>2000){
                //params.reduce_w();
                params.relax_tolerance();
                not_cnv_cnt = 0;
            }

            //Linear mixing for n and p
            for(int i = 1;i<num_cell;i++){
                p[i] = newp[i]*params.w + oldp[i]*(1.-params.w);
                n[i] = newn[i]*params.w + oldn[i]*(1.-params.w);
            }
            iter = iter+1;
        }

        std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(finish-start);
        std::cout << "1 Va CPU time = " << time.count() << std::endl;

        //-------------------Calculate Currents using Scharfetter-Gummel definition--------------------------

        p[0] = continuity_p.get_p_leftBC();
        n[0]  = continuity_n.get_n_leftBC();
        for(int i = 1; i<num_cell;i++){
            Jp[i] = -(q*Vt*params.N*params.mobil/params.dx)*p_mob[i]*(p[i]*B_p2[i]-p[i-1]*B_p1[i]);
            Jn[i] =  (q*Vt*params.N*params.mobil/params.dx)*n_mob[i]*(n[i]*B_n1[i]-n[i-1]*B_n2[i]);
            J_total[i] = Jp[i] + Jn[i];
        }

        //---------------------Write current Va step's data to file---------------------------------------------

        //Write charge densities, recombination rates, etc
        std::string filename = std::to_string(Va);
        filename += ".txt";  //add .txt extension
        VaData.open(filename); //this will need to have a string as file name
        for(int i = 1;i<=num_cell;i++){
            VaData << std::setw(15) << std::setprecision(8) << params.dx*i;
            VaData << std::setw(15) << std::setprecision(8) << Vt*V[i];
            VaData << std::setw(15) << std::setprecision(8) << params.N*p[i];         //setprecision(8) sets that use 8 sigfigs
            VaData << std::setw(15) << std::setprecision(8) << params.N*n[i];
            VaData << std::setw(15) << std::setprecision(8) << J_total[i];
            VaData << std::setw(15) << std::setprecision(8) << Un[i];
            VaData << std::setw(15) << std::setprecision(8) << PhotogenRate[i];
            VaData << std::setw(15) << std::setprecision(8) << R_Langevin[i];
            VaData << std::setw(15) << std::setprecision(8) << params.w;
            VaData << std::setw(15) << std::setprecision(8) << params.tolerance <<std::endl;
        }
        VaData.close();

        //Write JV file
        if(Va_cnt >0){
            if(JV.is_open()) {
                JV << Va << " " << J_total[num_cell-1] << " " << iter << "\n";
            }
        }
    }
    JV.close();

    return 0;
}
