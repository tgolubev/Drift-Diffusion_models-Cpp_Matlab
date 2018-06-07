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
%     for photogeneration rate can be added to photogeneration.cpp. Generation rate file
%     should contain num_cell -2 number of entries in a single column, corresponding to
%     the the generation rate at each mesh point (except the endpoints).
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
#include "simulation.h"
#include "poisson.h"
#include "continuity_p.h"
#include "continuity_n.h"
#include "recombination.h"
#include "photogeneration.h"
#include "thomas_tridiag_solve.h"
#include "bernoulli.h"


int main()
{
    std::ofstream JV;
    std::ofstream VaData;
    JV.open("JV.txt");  //note: file will be created inside the build directory

    Simulation simul;   //constuct Simulation object (this will initialize all simulation parameters, in future will read from file)

    double old_error;

    //Fill the RELATIVE dielectric constant vector (can be position dependent, as long as piecewise constant)
    std::vector<double> epsilon(simul.get_num_cell()+1);
    std::fill(epsilon.begin(), epsilon.end(), eps_active);

    //Fill the mobilities vectors (can be position dependent, as long as piecewise constant)
    std::vector<double> p_mob(simul.get_num_cell()+1), n_mob(simul.get_num_cell()+1);
    std::fill(p_mob.begin(), p_mob.end(), p_mob_active);
    std::fill(n_mob.begin(), n_mob.end(), n_mob_active);

    for(int i = 0;i<= simul.get_num_cell();i++){
        p_mob[i] = p_mob[i]/mobil;
        n_mob[i] = n_mob[i]/mobil;
    }

    //Initialize other vectors
    //Note: vectors are automatically initialized to 0.
    //Will use indicies for n and p... starting from 1 --> since is more natural--> corresponds to 1st node inside the device...
    std::vector<double> n(simul.get_num_cell()), p(simul.get_num_cell()), oldp(simul.get_num_cell()), newp(simul.get_num_cell()), oldn(simul.get_num_cell()), newn(simul.get_num_cell()), oldV(simul.get_num_cell()+1), newV(simul.get_num_cell()+1), V(simul.get_num_cell()+1);
    std::vector<double> Un(simul.get_num_cell()), Up(simul.get_num_cell()), R_Langevin(simul.get_num_cell()), PhotogenRate(simul.get_num_cell());  //store the results of these..
    std::vector<double> Jp(simul.get_num_cell()),Jn(simul.get_num_cell()), J_total(simul.get_num_cell());

    //DO SOMETHING WITH BERNOULLIS--> PUT  AS PART OF CONTINUITY CLASSES?
    std::vector<double> B_n1(simul.get_num_cell()+1), B_n2(simul.get_num_cell()+1), B_p1(simul.get_num_cell()+1), B_p2(simul.get_num_cell()+1); //vectors for storing Bernoulli fnc values

    std::vector<double> error_np_vector(simul.get_num_cell());  //for storing errors between iterations

    //Boundary conditions (FOR NOW OK, READS FROM parameters.h)
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
    double diff = (V_rightBC-V_leftBC)/simul.get_num_cell();
    V[0] = V_leftBC;  //fill V(0) here for use in Beroulli later
    for(int i = 1; i<=simul.get_num_cell()-1; i++){
        V[i] = V[i-1] + diff;
    }
    V[simul.get_num_cell()] = V_rightBC;


    //////////////////////MAIN LOOP////////////////////////////////////////////////////////////////////////////////////////////////////////
    int iter;
    double error_np;
    bool not_converged;
    int not_cnv_cnt;
    int Va_cnt;
    double Va;


    Poisson poisson(simul);  //set up a poisson eqn object
    Recombo recombo(simul, k_rec, E_trap, n1, p1);  //set up recombination object
    Continuity_n continuity_n(simul);
    Continuity_p continuity_p(simul);
    Photogeneration photogen(simul, Photogen_scaling);

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();  //start clock timer
    for(Va_cnt = 0; Va_cnt <=simul.get_num_V() +1;Va_cnt++){  //+1 b/c 1st Va is the equil run
        not_converged = false;
        not_cnv_cnt = 0;
        if(simul.get_tolerance() > 1e-5){
            std::cout<<"ERROR: Tolerance has been increased to > 1e-5" <<std::endl;
        }
        if(Va_cnt==0){
            simul.use_tolerance_eq();  //relaxed tolerance for equil. run
            simul.use_w_eq();
            Va = 0;
        }
        else{
            Va = simul.get_Va_min()+simul.get_increment()*(Va_cnt-1);
        }
        if(Va_cnt ==1){
            simul.use_tolerance_i();  //reset tolerance back
            simul.use_w_i();
            PhotogenRate = photogen.getPhotogenRate();    //otherwise PhotogenRate is pre-initialized to 0 in this main.cpp when declared
        }
        std::cout << "Va = " << Va <<std::endl;

        //Apply the voltage boundary conditions
        V_leftBC = -((Vbi  -Va)/(2*Vt)-phi_a/Vt);
        V_rightBC = (Vbi- Va)/(2*Vt) - phi_c/Vt;
        //adjust them within V
        V[0] = V_leftBC;
        V[simul.get_num_cell()] = V_rightBC;

        error_np = 1.0;
        iter = 0;
        while(error_np > simul.get_tolerance()){
            std::cout << "error np " << error_np <<std::endl;
            std::cout << "Va " << Va <<std::endl;
            //-----------------Solve Poisson Equation------------------------------------------------------------------
            poisson.set_rhs(epsilon, n, p, V_leftBC, V_rightBC); //setup RHS of Poisson eqn (bV)
            oldV = V;
            newV = Thomas_solve(poisson.get_main_diag(), poisson.get_upper_diag(), poisson.get_lower_diag(), poisson.get_rhs());

            //add on the BC's --> b/c matrix solver just outputs the insides...
            newV[0] = V[0];
            newV[simul.get_num_cell()] = V[simul.get_num_cell()];

            //Mix old and new solutions for V
            if(iter>0){
                for(int i =  1;i<=simul.get_num_cell()-1;i++){
                    V[i] = newV[i]*simul.get_w() + oldV[i]*(1.-simul.get_w());
                }
            }
            else V = newV;

            //------------------------------Calculate Net Generation Rate----------------------------------------------------------
            R_Langevin = recombo.ComputeR_Langevin(simul, n,p);

            for(int i = 1;i<=simul.get_num_cell()-1;i++){
                Un[i] = PhotogenRate[i] - R_Langevin[i];
            }
            Up = Un;

            //--------------------------------Solve equations for n and p------------------------------------------------------------ 
            //setup the matrices
            BernoulliFnc_n(simul, V, B_n1, B_n2);  //fills B_p1 and B_p2 with the values of Bernoulli fnc
            continuity_n.set_main_diag(n_mob, B_n1, B_n2);
            continuity_n.set_upper_diag(n_mob, B_n1);
            continuity_n.set_lower_diag(n_mob, B_n2);
            continuity_n.set_rhs(n_mob, B_n1, B_n2, n_leftBC, n_rightBC, Un);
            oldn = n;
            newn = Thomas_solve(continuity_n.get_main_diag(), continuity_n.get_upper_diag(), continuity_n.get_lower_diag(), continuity_n.get_rhs());

            BernoulliFnc_p(simul, V, B_p1, B_p2);
            continuity_p.set_main_diag(p_mob, B_p1, B_p2);
            continuity_p.set_upper_diag(p_mob, B_p2);
            continuity_p.set_lower_diag(p_mob, B_p1);
            continuity_p.set_rhs(p_mob, B_p1, B_p2, p_leftBC, p_rightBC, Up);
            oldp = p;
            newp = Thomas_solve(continuity_p.get_main_diag(), continuity_p.get_upper_diag(), continuity_p.get_lower_diag(), continuity_p.get_rhs());

            //if get negative p's or n's set them = 0
            for(int i = 1; i<=simul.get_num_cell()-1;i++){
                if(newp[i]<0.0) newp[i] = 0;
                if(newn[i]<0.0) newn[i] = 0;
            }

            //calculate the error
            old_error = error_np;
            for (int i = 1;i<=simul.get_num_cell()-1;i++){
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
                simul.relax_tolerance();
                not_cnv_cnt = 0;
            }

            //Linear mixing for n and p
            for(int i = 1;i<=simul.get_num_cell()-1;i++){
                p[i] = newp[i]*simul.get_w() + oldp[i]*(1.-simul.get_w());
                n[i] = newn[i]*simul.get_w() + oldn[i]*(1.-simul.get_w());
            }
            iter = iter+1;
        }

        std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(finish-start);
        std::cout << "1 Va CPU time = " << time.count() << std::endl;

        //-------------------Calculate Currents using Scharfetter-Gummel definition--------------------------
        p[0] = p_leftBC;
        n[0]  = n_leftBC;
        for(int i = 1; i<=simul.get_num_cell()-1;i++){
            Jp[i] = -(q*Vt*N*mobil/dx)*p_mob[i]*(p[i]*B_p2[i]-p[i-1]*B_p1[i]);
            Jn[i] =  (q*Vt*N*mobil/dx)*n_mob[i]*(n[i]*B_n1[i]-n[i-1]*B_n2[i]);
            J_total[i] = Jp[i] + Jn[i];
        }

        //---------------------Write current Va step's data to file---------------------------------------------
        //Write charge densities, recombination rates, etc
        std::string filename = std::to_string(Va);
        filename += ".txt";  //add .txt extension
        VaData.open(filename); //this will need to have a string as file name
        for(int i = 1;i<=simul.get_num_cell();i++){
            VaData << std::setw(15) << std::setprecision(8) << dx*i;
            VaData << std::setw(15) << std::setprecision(8) << Vt*V[i];
            VaData << std::setw(15) << std::setprecision(8) << N*p[i];         //setprecision(8) sets that use 8 sigfigs
            VaData << std::setw(15) << std::setprecision(8) << N*n[i];
            VaData << std::setw(15) << std::setprecision(8) << J_total[i];
            VaData << std::setw(15) << std::setprecision(8) << Un[i];
            VaData << std::setw(15) << std::setprecision(8) << PhotogenRate[i];
            VaData << std::setw(15) << std::setprecision(8) << R_Langevin[i];
            VaData << std::setw(15) << std::setprecision(8) << simul.get_w();
            VaData << std::setw(15) << std::setprecision(8) << simul.get_tolerance() <<std::endl;
        }
        VaData.close();

        //Write JV file
        if(Va_cnt >0){
            if(JV.is_open()) {
                JV << Va << " " << J_total[simul.get_num_cell()-1] << " " << iter << "\n";  //NOTE: still some issue here, floor(num_cell/2), in the middle get different currents
            }
        }
    }
    JV.close();



    return 0;
}
