/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving 1D Poisson + Drift Diffusion semiconductor eqns using
%                    Scharfetter-Gummel discretization
%
%                  Coded by Timofey Golubev (2018.05.26)
%               NOTE: i=1 corresponds to x=0, i=nx to x=L
%
%     The code as is will calculate and plot a JV curve
%     as well as carrier densities, current densities, and electric field
%     distributions of a generic solar cell. More equations for carrier
%     recombination can be added. Generation rate will be inputted from gen_rate.txt
%     file (i.e. the output of an optical model can be used) or an analytic expression
%     for photogeneration rate can be added to photogeneration.cpp.
%
%     The code can also be applied to any semiconductor device by
%     setting photogeneration rate to 0 and adding equations for
%     loss mechanisms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>   //allows to use fill and min
#include <fstream>
#include <chrono>
#include <string>
#include <time.h>
#include "parameters.h"
#include "Set_diagonals.h"
#include "bernoulli.h"
#include "recombination.h"
#include "set_rhs.h"
#include "photogeneration.h"
#include "thomas_tridiag_solve.h"

int main()
{
    std::ofstream JV;
    std::ofstream VaData;
    JV.open("JV.txt");  //note: file will be created inside the build directory
    double Va;
    double w;
    double tolerance;
    double old_error;

    //Fill the RELATIVE dielectric constant vector (can be position dependent, as long as piecewise constant)
    std::vector<double> epsilon(num_cell+1);
    std::fill(epsilon.begin(), epsilon.end(), eps_active);

    //Fill the mobilities vectors (can be position dependent, as long as piecewise constant)
    std::vector<double> p_mob(num_cell+1), n_mob(num_cell+1);
    std::fill(p_mob.begin(), p_mob.end(), p_mob_active);
    std::fill(n_mob.begin(), n_mob.end(), n_mob_active);

    for(int i = 0;i<= num_cell;i++){
        p_mob[i] = p_mob[i]/mobil;
        n_mob[i] = n_mob[i]/mobil;
    }

    //Initialize other vectors
    //Note: vectors are automatically initialized to 0.
    //Will use indicies for n and p... starting from 1 --> since is more natural--> corresponds to 1st node inside the device...
    std::vector<double> n(num_cell), p(num_cell), oldp(num_cell), newp(num_cell), oldn(num_cell), newn(num_cell), oldV(num_cell+1), newV(num_cell+1), V(num_cell+1);
    std::vector<double> Un(num_cell), Up(num_cell), Photogen_rate(num_cell); //bn, bp are rhs of continuity eqns
    std::vector<double> a (num_cell);//main diag
    std::vector<double> b(num_cell-1); //upper diag, size can be = num_elements b/c off-diags are 1 element less than  main diag
    std::vector<double> c(num_cell-1);//lower diag
    std::vector<double> a_Poisson(num_cell), b_Poisson(num_cell-1), c_Poisson(num_cell-1); //for Poisson matrix diagonals
    std::vector<double> rhs(num_cell);
    std::vector<double> Jp(num_cell),Jn(num_cell), J_total(num_cell);
    std::vector<double> B_n1(num_cell+1), B_n2(num_cell+1), B_p1(num_cell+1), B_p2(num_cell+1); //vectors for storing Bernoulli fnc values
    std::vector<double> R_Langevin(num_cell), R_SRH_HTL(num_cell), R_SRH_ETL(num_cell);
    std::vector<double> error_np_vector(num_cell);  //for storing errors between iterations

    //Boundary conditions
    double n_leftBC = (N_LUMO*exp(-(E_gap - phi_a)/Vt))/N;       //this is anode
    double p_leftBC = (N_HOMO*exp(-phi_a/Vt))/N;
    double n_rightBC = (N_LUMO*exp(-phi_c/Vt))/N;
    double p_rightBC = (N_HOMO*exp(-(E_gap - phi_c)/Vt))/N;

    //Initial conditions
    double min_dense = std::min(n_leftBC,  p_rightBC);
    std::fill(n.begin()+1, n.end(), min_dense);
    std::fill(p.begin()+1, p.end(), min_dense);

    double V_leftBC = -((Vbi)/(2*Vt)-phi_a/Vt);
    double V_rightBC = (Vbi)/(2*Vt)-phi_c/Vt;
    double diff = (V_rightBC-V_leftBC)/num_cell;
    V[0] = V_leftBC;  //fill V(0) here for use in Beroulli later
    for(int i = 1; i<=num_cell-1; i++){
        V[i] = V[i-1] + diff;
    }
    V[num_cell] = V_rightBC;

    //////////////////////MAIN LOOP////////////////////////////////////////////////////////////////////////////////////////////////////////
    int iter;
    double error_np;
    bool not_converged;
    int not_cnv_cnt;
    int Va_cnt;

    //Setup Poisson matrix (need to do only once)
    set_main_AVdiag(epsilon, a_Poisson);
    set_upper_AVdiag(epsilon, b_Poisson);
    set_lower_AVdiag(epsilon, c_Poisson);

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();  //start clock timer
    for(Va_cnt = 0; Va_cnt <=num_V +1;Va_cnt++){  //+1 b/c 1st Va is the equil run
        not_converged = false;
        not_cnv_cnt = 0;
        if(tolerance > 1e-5){
            std::cout<<"ERROR: Tolerance has been increased to > 1e-5" <<std::endl;
        }
        if(Va_cnt==0){
            tolerance = tolerance_eq;  //relaxed tolerance for equil. run
            w = w_eq;
            Va = 0;
        }
        else{
            Va = Va_min+increment*(Va_cnt-1);
        }
        if(Va_cnt ==1){
            tolerance =  tolerance_i;  //reset tolerance back
            w = w_i;
            Photogen_rate = PhotogenerationRate();  //ONLY NEED TO DO THIS ONCE FOR 1ST nonequil Va, since photogen rate stays constant!
        }
        std::cout << "Va = " << Va <<std::endl;

        //Apply the voltage boundary conditions
        V_leftBC = -((Vbi  -Va)/(2*Vt)-phi_a/Vt);
        V_rightBC = (Vbi- Va)/(2*Vt) - phi_c/Vt;
        //adjust them within V
        V[0] = V_leftBC;
        V[num_cell] = V_rightBC;

        error_np = 1.0;
        iter = 0;
        while(error_np > tolerance){
            std::cout << "error np " << error_np <<std::endl;
            std::cout << "Va " << Va <<std::endl;
            //-----------------Solve Poisson Equation------------------------------------------------------------------
            set_bV(epsilon, n, p, V_leftBC, V_rightBC, rhs); //setup RHS of Poisson eqn (bV)
            oldV = V;
            newV = Thomas_solve(a_Poisson, b_Poisson, c_Poisson, rhs);

            //add on the BC's --> b/c matrix solver just outputs the insides...
            newV[0] = V[0];
            newV[num_cell] = V[num_cell];

            //Mix old and new solutions for V
            if(iter>0){
                for(int i =  1;i<=num_cell-1;i++){
                    V[i] = newV[i]*w + oldV[i]*(1.-w);
                }
            }
            else V = newV;

            //------------------------------Calculate Net Generation Rate----------------------------------------------------------
            ComputeR_Langevin(n,p, R_Langevin);

            for(int i = 1;i<=num_cell-1;i++){
                Un[i] = Photogen_rate[i] - R_Langevin[i];
            }
            Up = Un;


            //--------------------------------Solve equations for n and p------------------------------------------------------------ 
            //setup the matrices
            BernoulliFnc_n(V, B_n1, B_n2);  //fills B_p1 and B_p2 with the values of Bernoulli fnc
            set_main_An_diag(n_mob, B_n1, B_n2, a);
            set_upper_An_diag(n_mob, B_n1, b);
            set_lower_An_diag(n_mob, B_n2, c);
            set_bn(n_mob, B_n1, B_n2, n_leftBC, n_rightBC, Un, rhs);
            oldn = n;
            newn = Thomas_solve(a, b, c, rhs);

            BernoulliFnc_p(V, B_p1, B_p2);
            set_main_Ap_diag(p_mob, B_p1, B_p2, a);
            set_upper_Ap_diag(p_mob, B_p2, b);
            set_lower_Ap_diag(p_mob, B_p1, c);
            set_bp(p_mob, B_p1, B_p2, p_leftBC, p_rightBC, Up, rhs);
            oldp = p;
            newp = Thomas_solve(a, b, c, rhs);

            //if get negative p's or n's set them = 0
            for(int i = 1; i<=num_cell-1;i++){
                if(newp[i]<0.0) newp[i] = 0;
                if(newn[i]<0.0) newn[i] = 0;
            }

            //calculate the error
            old_error = error_np;
            for (int i = 1;i<=num_cell-1;i++){
                if(newp[i]!= 0 && newn[i] !=0){
                    error_np_vector[i] = (abs(newp[i]-oldp[i])+abs(newn[i]-oldn[i]))/abs(oldp[i]+oldn[i]);
                }
            }
            error_np = *std::max_element(error_np_vector.begin(),error_np_vector.end());
            std::fill(error_np_vector.begin(), error_np_vector.end(),0.0);  //refill with 0's so have fresh one for next iter

            //auto decrease w if not converging
            if(error_np>=old_error) not_cnv_cnt = not_cnv_cnt+1;
            if(not_cnv_cnt>2000){
                //w = w/2.  //don't auto decrease for 0nm C60 case
                tolerance = tolerance*10;
                not_cnv_cnt = 0;
            }

            //Linear mixing for n and p
            for(int i = 1;i<=num_cell-1;i++){
                p[i] = newp[i]*w + oldp[i]*(1.-w);
                n[i] = newn[i]*w + oldn[i]*(1.-w);
            }
            iter = iter+1;
        }

        std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(finish-start);
        std::cout << "1 Va CPU time = " << time.count() << std::endl;

        //-------------------Calculate Currents using Scharfetter-Gummel definition--------------------------
        p[0] = p_leftBC;
        n[0]  = n_leftBC;
        for(int i = 1; i<=num_cell-1;i++){
            Jp[i] = -(q*Vt*N*mobil/dx)*p_mob[i]*(p[i]*B_p2[i]-p[i-1]*B_p1[i]);
            Jn[i] =  (q*Vt*N*mobil/dx)*n_mob[i]*(n[i]*B_n1[i]-n[i-1]*B_n2[i]);
            J_total[i] = Jp[i] + Jn[i];
        }

        //---------------------Write current Va step's data to file---------------------------------------------
        //Write charge densities, recombination rates, etc
        std::string filename = std::to_string(Va);
        VaData.open(filename); //this will need to have a string as file name
        for(int i = 1;i<num_cell-1;i++){
            VaData << std::setw(15) << std::setprecision(8) << dx*i;
            VaData << std::setw(15) << std::setprecision(8) << Vt*V[i];
            VaData << std::setw(15) << std::setprecision(8) << N*p[i];         //setprecision(8) sets that use 8 sigfigs
            VaData << std::setw(15) << std::setprecision(8) << N*n[i];
            VaData << std::setw(15) << std::setprecision(8) << J_total[i];
            VaData << std::setw(15) << std::setprecision(8) << Un[i];
            VaData << std::setw(15) << std::setprecision(8) << Photogen_rate[i];
            VaData << std::setw(15) << std::setprecision(8) << R_Langevin[i];
            VaData << std::setw(15) << std::setprecision(8) << R_SRH_ETL[i];
            VaData << std::setw(15) << std::setprecision(8) << R_SRH_HTL[i];
            VaData << std::setw(15) << std::setprecision(8) << w;
            VaData << std::setw(15) << std::setprecision(8) << tolerance <<std::endl;
        }
        VaData.close();

        //Write JV file
        if(Va_cnt >0){
            if(JV.is_open()) {
                JV << Va << " " << J_total[num_cell-1] << " " << iter << "\n";  //NOTE: still some issue here, floor(num_cell/2), in the middle get different currents
            }
        }
    }
    JV.close();

    return 0;
}
