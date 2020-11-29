%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving Poisson + Drift Diffusion eqns (for holes) using
%                   Scharfetter-Gummel discretization
%                        with Newton's method
%
%                  Coded by Timofey Golubev (2017.10.07)
%             NOTE: i=1 corresponds to x=0, i=nx to x=L

%RECALL THAT fullV = V/Vt: is scaled!. Same with p = p/N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;   %NOTE: clear improves performance over clear all, and still clears all variables.

%% Declare global variables  (accessible by all functions)
global L l num_cell num_elements Cp Cn  N q kb T Vt dx L_int epsilon p_mob n_mob CV CnT CpT Cp_gummel Cn_gummel mobil;
global H_A H_D T_tD T_tA T N_LUMO N_HOMO gap_A gap_D cap_n cap_p Eb rc k_rec phi_c phi_a gap_elec gap_hole;
global k_ppd_eq k_ppr_0 Jx %this is needed for RecalculateU function


%% Parameters
L = 73*10^-9;             %device length in meters
l = 51;                   %position of interface
num_cell = 73;            % number of cells

p_mob = 2.0*10^-8;        %hole mobility
n_mob =  2.0*10^-7;
mobil =  2.0*10^-7;       %scaling for mobility: NOTE: is used for Gummel but not for Newton

dx = L/num_cell;
L_int = dx;               %Interface width, NOTE: if make L_int!= dx, must change bV(l) and l+1 calculations to reflect this!

N_LUMO = 6.0*10^27;
N_HOMO = 6.0*10^27;

H_A = 1.0D24;			  %Trap DOS in acceptor (1/m^3)
H_D = 1.0D24;
%Energy levels (Currently for P3HT:PCBM)
LUMO_D = -3.0;
HOMO_D = -4.9;
LUMO_A = -3.7;
HOMO_A = -6.1;
phi_a = 0.2;			%Anode injection barrier
phi_c = 0.1;		    %Cathode injection barrier

Jx = 4.4D20;       %MAKING THIS HIGHER, JUST MAKES u BECOME NEGATIVE AT SLIGHTLY HIGHER VOLTAGE U %exciton current density --> for calculating zeta : PP density\

%calculated parameters
E_gap = LUMO_A - HOMO_D;        %Interface energy gap
gap_D = LUMO_D - HOMO_D;		%Donor band gap
gap_A = LUMO_A - HOMO_A;		%Acceptor band gap
gap_hole = HOMO_D - HOMO_A;		%Interface hole energy offset
gap_elec = LUMO_D - LUMO_A;		%Interface electron energy offset

Vbi = E_gap - phi_a - phi_c;      %Built-in voltage: this controlls when current in solar cell will switch direction


N = 10^27.;                    %scaling factor helps CV be on order of 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Va_min = 0.00;               %volts
Va_max = 1.2;
increment = 0.01;            %for increasing V
num_V = floor((Va_max-Va_min)/increment)+1;
%num_V = 200;

%Simulation parameters
w_i = 0.2;                    %set up of weighting factor: for 1st Va--> Gummel is used
w_newton = 0.6;  %NOTE: IF GETTING 0'S in the charge densities, means w is too big--> have oscillating soln's giving negative densities which I set to 0!
%and higher w for rest.
tolerance = 5*10^-13;        %error tolerance
tolerance_i = 5*10^-13;
equil_tolerance = 5*10^-12;  %is for Gummel --> equil convergence
manual_Rp = true;
negative_nps = false; %for checks of whether have negative n's and/or p's after update in Newton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Physical Constants
q =  1.60217646*10^-19;         %elementary charge, C
kb = 1.3806503*10^-23;          %Boltzmann const., J/k
T = 296.;                       %temperature
T_tD = 680.0D0;			        %Trap temperature (K)
T_tA = 680.0D0;
epsilon_0 =  8.85418782*10^-12; %F/m
epsilon = 3.8*epsilon_0;        %dielectric constant of P3HT:PCBM

Vt = (kb*T)/q;

%Recombination parameters
k_ppr_0 = (6.0D-5)^(-1.0D0);                %NOTE: matlab understands the D notation
rc = q^2/(4.0D0*pi*epsilon*kb*T);	        %Exiton separation length (Coloumb radius)
Eb  = q/(4.0D0*pi*epsilon*L_int);   	    %Exciton binding energy
Rp_man = 2000.;                             %manual parallel resistence

k_rec = q*(n_mob+ p_mob)/epsilon;      %from Langevin
k_rec_n = q*n_mob/epsilon;
k_rec_p = q*p_mob/epsilon;

cap_n = k_rec_n;   %these are used in f_Dr_SRH functions
cap_p = k_rec_p;

%% Domain Discretization
a_1=0; a_2=L; x=linspace(a_1,a_2,num_cell+1); %dx=(b-a)/num_cell;   %x is positions array, +1 necessary b/c linspace needs (# of pts = # of cells +1)
nx = length(x);
num_elements = nx-2;

%% Boundary Conditions
n_full(1) = N_LUMO*exp(-gap_D/Vt);
n_full(num_cell+1) = N_LUMO;
p_full(1) = N_HOMO;
p_full(num_cell+1) = N_HOMO*exp(-gap_A/Vt);

n_full = n_full/N;
p_full = p_full/N;

%% Initial Conditions

min_dense = min(n_full(1),  p_full(num_cell+1)); %couldn't figure out how to take min  of all 4 bc's so just take min of the 2 bc's that I know are lowest --> opposite side of their electrodes.

for i = 1:num_elements         %initial conditions for QR etc. %NOTE: p(1) = p_full(2)  etc.
    p(i) = min_dense;
    n(i) = min_dense;
end

%set up initial guess of V for QR algorithm:
%make V linearly decreasing
fullV(1)          = 0.0;
fullV(num_cell+1) = (Vbi)/Vt;  %for intial guess: Va = 0, b/c is equil. run, so Vbi-Va --> Vbi
diff = (fullV(num_cell+1)-fullV(1))/num_cell;   %note I have 1 less pts than in the fortran version
for i = 2:num_elements+1
    fullV(i) = fullV(i-1)+diff;
end
V_initial = fullV(2:num_elements+1);  %shift to correspond to V in  matrix: this will be used as initial guess in QR solver
V_initial = V_initial.';   %transpose to column vector

%% Main voltage loop
Va_cnt = 0;
clear p_solution;
clear iterations;
for Va_cnt = 0:num_V+1   %+1 b/c 1st Va is the equil. Va = 0 run
    not_converged = false;  %everytime gets to here, means it has converged, so make this false, so tolerance not relaxed unless not converged occurs at least once
    
    %stop the calculation if tolerance becomes too high
    if(tolerance >10^-4)
        break
    end
    
    not_cnv_cnt = 0;
    %allocate matrices/arrays
    V = zeros(nx-2,1);  %note: allocate this before set V = V_initial for eq run!
    
    if(Va_cnt ==0)
        tolerance = equil_tolerance;       %relax tolerance for 1st convergence
        w = 0.01;   %set w to be small for 1st convergence: like in fortran version
        Va = 0;
        Un = 0;   %U is always for 1st Va --> this is equilibrium calculation
        Up = 0;
        V = V_initial;
    end
    if(Va_cnt ==1)
        tolerance = tolerance_i;       %reset tolerance back
        w = w_i;
    end
    
    if(Va_cnt >0 )
        Va = Va_min+increment*(Va_cnt-1);     %set Va value
        %Va = Va_max-increment*(Va_cnt-1);    %decrease Va by increment in each iteration
    end
    
    fullV(1)          = 0.0;
    fullV(num_cell+1) = (Vbi -Va)/Vt;
    
    %Need to define these here b/c need to use p(l) in Rpmanual below!:
    %otherwise in 1st run, will not have a value for p(l) etc...
    p_full = [p_full(1), p, p_full(num_cell+1)]; %bndrys already defined above: need ; FOR IF NOT USING QR FOR BERNOULLI EQNS
    n_full = [n_full(1), n, n_full(num_cell+1)];
    
    %timing
    tic
    
    %% Setup the solver
    iter = 1;
    error_np =  1.0;
    Jp = zeros(1,num_cell);
    Jn = zeros(1,num_cell);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solver Loop
    %NOTE: Cp and Cn refer to the coeffs for Newtons iterations...
    Cp = dx^2/(Vt*N*p_mob);   %note:now use mobil: the scaling here: b/c the actual p and n mob are inside the matrix
    Cn = dx^2/(Vt*N*n_mob);
    Cp_gummel = dx^2/(Vt*N*mobil);
    Cn_gummel = dx^2/(Vt*N*mobil);
    CV = N*dx^2*q/(epsilon*Vt);
    CnT = (H_A/N_LUMO)*(T/T_tA);  %coefficient in front of dFnT/dn derivative
    CpT = (H_D/N_HOMO)*(T/T_tD);
    
    %I NEED TO START THE SOLVER WITH GUMMEL METHOD!! --> TO GET INITIAL
    %GUESS TO BE CLOSE, AND THEN CONTINUE WITH NEWTON!!
    
    %%  Gummel solver (to begin with)
    %set up AV: these never change so do outside of loop
    %set up using sparse matrix (just store the lower diag, main diag,upper
    %diag in num_pts x 3 matrix AV_val.
    
    %Here don't have to worry about the exactl setup of colunms for sparse
    %matrix b/c all values in each diagonal are same anyway.
    AV = zeros(nx-2,3);
    AV(:,1) = 1.;
    AV(:,3) = 1.;
    AV(:,2) = -2.;
    AV_val = spdiags(AV,-1:1,nx-2,nx-2);  %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
    bV = zeros(nx-2,1);
    
    if(Va_cnt == 0)
        num_Gummel_iters = 2800;
    else
        num_Gummel_iters =  120;
    end
    
    
    %      while (error_np > tolerance && iter < num_Gummel_iters )
    %base exit out of loop on how large the error is--> MAKE SURE ERROR IS
    %SMALL ENOUGH BEFORE ENTER Newton loop...
    
    %ONLY USE GUMMEL  ON equil iter and 1st non-equil Va (b/c there we add
    %a Un and Up
    if(Va_cnt <2)
        
        while (error_np > tolerance)  %takes gummel a little more than 2900 iters to converge, so want to switch to Newton just prior to that for testing
            
            %WE WILL call traps here, but not in Newton!--> there just use
            %the traps from last Gummel iter--> traps don't change much
            %anyway!
            
            [pT_full, nT_full] = CalTrap(p_full, n_full);   %calculate trap densities
            nT = nT_full(2:num_cell);   %redefine as only the points inside device
            pT = pT_full(2:num_cell);
            
            %% Poisson equation with tridiagonal solver
            % setup bV
            for i = 1: num_elements
                bV(i,1) = CV*((n(i)+nT(i))-(p(i)+pT(i)));   %THIS IS THE VERSION WITHOUT TRAPS--> TO SIMPLIFY THE JACOBIAN
            end
            %for bndrys
            bV(1,1) = bV(1,1) - fullV(1);
            bV(num_elements,1) = bV(num_elements,1) - fullV(num_cell+1);
            
            %call solver, solve for V
            %spparms('spumoni',2)     %THIS ALLOWS TO SEE THE DETAILS OF THE
            %SOLVER CONDITIONS!
            
            oldV = V;
            newV = lsqr(AV_val,bV,10^-15,1000,[] ,[] ,V);  %can't converge to below 10^-14 on the QR here
            
            normV = norm(AV_val*V-bV)  %check convergence again
            if(iter >1)
                V = newV*w + oldV*(1.-w);
            else
                V = newV;  %no mixing for 1st iter
            end
            
            %V =  AV_val\bV;
            %SO THIS ACTUALLY IS SOLVING FOR PSI PRIME: A SCALED PSI! --> so
            %later in Bernoulli: don't need to devide by VT again: b/c this
            %here is already scaled!
            
            %make V with bndry pts
            fullV = [fullV(1); V; fullV(num_cell+1)];       %note: bndry values are already defined above
            
            fullV = fullV.';             %transpose
            
            %NOTE: IT IS WRONG TO SCALE fullV so don't do that!
            %------------------------------------------------------------------------------------------------
            %% Calculate recombination rates
            % Calculate recombination rates pn each iteration within the Gummel
            % steps -->  b/c works, but in Newton, don't calculate them at
            % all!! --> should be already pre-converged!! from Gummel
            if(Va_cnt >0)  %for all but the first Va, we calculate the rates: 1st Va is equil. run
                F = (fullV(l+1)-fullV(l)).*Vt./L_int;
                k_ppd =  kppd(F);
                Rtotal = R(n_full, p_full);  %calculates recombination:  we use trap and SRH so implement just those 2--> do in  sepearte function file
                k_ppr = k_ppr_0.*exp(-F.*L_int./Vt);
                %test = f_of_zeta(pi/2, F)   %f of zeta function works
                zeta_eq = Rtotal./k_ppd_eq;   %use ./ to tell it that these are not matrixes! otherwise get error.
                zeta = (Jx./L_int+ k_ppr.*zeta_eq + Rtotal)./(k_ppr+k_ppd);
                
                Un = k_ppd.*zeta - Rtotal;
                Up = k_ppd.*zeta - Rtotal;
                if(manual_Rp)
                    temp = ( E_gap - F*L_int+ Vt*log(n_full(l+1)*p_full(l)*N^2/(N_LUMO*N_HOMO)) )/(q*L_int*Rp_man);  %IN TEH FIRST TRIAL, N_FULL L+1 IS NOT DEFINED!!! IS JUST = 0!!: NEED TO FIX THIS!!
                    Un = (Un - temp);   %temp is the manually inserted Rp value
                    Up = (Up - temp);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% now solve eqn for n and p
            
            %IN GUMMEL iters will use scaled mobilities just like before
            %scale mobilities
            p_mob = p_mob/mobil;
            n_mob = n_mob/mobil;
            
            B_n = BernoulliFnc_n(fullV);
            B_p = BernoulliFnc_p(fullV);
            Ap_val = SetAp_val(B_p);
            An_val = SetAn_val(B_n);
            bp = Setbp(B_p, p_full, Up);
            bn = Setbn(B_n, n_full, Un);
            
            old_p = p;
            p_sol = Ap_val\bp;
            normp = norm(Ap_val*p_sol-bp)
            % p_sol = lsqr(Ap_val,bp,10^-14,1000,[] ,[] ,p);
            %newp = p_sol;    %don't need to transpose b/c lsqr already gives column vector
            newp = p_sol.';  %need to transpose if don't use qr
            
            old_n = n;
            
            % n_sol = lsqr(An_val,bp,10^-14,1000,[] ,[] ,n);
            n_sol = An_val\bn;
            normn = norm(An_val*n_sol-bn)
            %newn = n_sol;
            newn = n_sol.';
            
            % if get negative p's or n's, make them equal 0
            %THIS IS FOR GUMMEL and is fine!
            for i = 1:num_elements
                if(newp(i) <0.0)
                    newp(i) = 0;
                end
                if(newn(i) <0.0)
                    newn(i) = 0;
                end
            end
            
            old_error =  error_np;
            count = 0;
            error_np_matrix = zeros(1,num_elements); %need to reset the matrix! b/c when more newp's become 0, need to remove that error element from matrix.
            for i = 1:num_elements
                if(newp(i) ~=0 && newn(i) ~=0)
                    count = count+1;  %counts number of non zero error calculations
                    error_np_matrix(count) = (abs(newp(i)-old_p(i))+abs(newn(i)-old_n(i)))/abs(old_p(i)+old_n(i));  %need the dot slash, otherwise it tries to do matrix operation! %ERROR SHOULD BE CALCULATED BEFORE WEIGHTING
                end
            end
            error_np = max(error_np_matrix)
            
            %auto decrease w if not converging
            if(error_np>= old_error)
                not_cnv_cnt = not_cnv_cnt+1;
            end
            if(not_cnv_cnt>1000)
                w = w/2.;
                %tolerance = tolerance*10;  %NOTE: this is Gummel tolerance--> which gets reduced if not converged-- only used at 1st iter--> Newton tolerance doesn't get reduced.
                not_cnv_cnt = 0;  %reset the count
            end
            w
            tolerance
            
            %weighting
            p = newp*w + old_p*(1.-w);
            n = newn*w + old_n*(1.-w);
            
            iter =  iter+1
            
            %unscale mobilities: using scaled mobilities for Gummel iters...
            p_mob = p_mob*mobil;
            n_mob = n_mob*mobil;
            
            
            %THESE NEED TO BE DEFINED INSIDE THE LOOP B/C I TAKE VALUES FROM PFULL
            %NFULL!!!
            p_full = [p_full(1), p, p_full(num_cell+1)]; %bndrys already defined above: need ; FOR IF NOT USING QR FOR BERNOULLI EQNS
            n_full = [n_full(1), n, n_full(num_cell+1)];
            
        end
        %save Un  and Up values so they can be passed to Newton solver:
        Un_gummel = Un;
        Up_gummel = Up;
    end
    
    
    % Calculate equil kppd so have it for future Va's
    if(Va_cnt == 0)
        F = (fullV(l+1)-fullV(l)).*Vt./L_int;
        k_ppd_eq = kppd(F);   %find equil. kppd value
    end
    
    
    %% Newton solver (will  NOT be used for equil run)
    %seems need to define V, since otherwise it equals 0!!
    V = fullV(2:num_cell);
    
    Newton_iter = 0;
    
    %file for outputing update_vector element for convergence behavior testing
    str = sprintf('%.3f',Va);
    filename = [str 'update_vector.txt']
    fupdate = fopen(fullfile('C:\Users\Tim\Documents\Duxbury_group_research\Drift_Diffusion_OLED\Newtons_method_OLED\results_output',filename),'w');   %w: Open or create new file for writing
    
    
    %make stricter tolerance for newton, for testing
    
    Newton_tolerance = 5*10^-13;  %can't converge to e-14 on higher Va's...
    
    
    %NOTE: to  not run  into issues related with calculating recombo rates,
    %which shouldn't be calculated in equil case, I CAN'T USE NEWTON for
    %equil calculation....
    
    if(Va_cnt >0)
        %reset the error value--> so start new with newton solver...
        error_np  = 1.0;
        while error_np > Newton_tolerance
            
            %timing
            tic
            
            Newton_iter = Newton_iter +1;
            if(Newton_iter == 1)
                
                %only call traps at 1st iteration!!.., so don't destabilize algo
                %too much
                [pT_full, nT_full] = CalTrap(p_full, n_full);   %calculate trap densities
                nT = nT_full(2:num_cell);   %redefine as only the points inside device
                pT = pT_full(2:num_cell);
                %Need to initialize these here, b/c Newton uses arrays for Un Up while
                %Gummel uses single values.
                Un = zeros(1, num_elements);
                Up = Un;
                %setting Un equal to Un from gummel is not needed b/c U's
                %are calculated below
                %                 Un(l) = Un_gummel;
                %                 Up(l-1) = Up_gummel;
                %redefine a general U, combining Un and Up:
                %                  U =  Un +  Up;
                
                %define solution column vector: will use the pre-partially converged
                %values from above Gummel method
                solution = [V.'; n.'; p.'];  %the transposes are so all are column vectors --> to correspond with update_vector...
                
            end
            
            B_n = BernoulliFnc_n(fullV);
            B_p = BernoulliFnc_p(fullV);
            
            % Calculate drift diffusion J's
            % Use the SG definition--> BUT PULL OUT THE COEFFS (and q's
            % cancel), and coeff's go to Cn and Cp
            for i = 1:num_cell        %verified that this is right
                %NOTE: NO coefficients, mobilities here, b/c those were moved
                %over to the Cn --> these J's are passed to Continuity  fnc's:
                %which use Jr - Jl +CnUn = 0 etc..., so coeffs are in the Cn
                Jp(1,i+1) = -(p_full(i+1)*B_p(2,i+1)-p_full(i)*B_p(1,i+1));      %(1,i)  b/c want a row vector  %need an N b/c my p's are scaled by N. Extra +1's are b/c B starts at i=2
                Jn(1,i+1) =  (n_full(i+1)*B_n(1,i+1)-n_full(i)*B_n(2,i+1));
            end
            
            
            if(Va_cnt >0)  %for all but the first Va, we calculate the rates: 1st Va is equil. run
                F = (fullV(l+1)-fullV(l)).*Vt./L_int;
                k_ppd =  kppd(F);
                Rtotal = R(n_full, p_full);  %calculates recombination:  we use trap and SRH so implement just those 2--> do in  sepearte function file
                k_ppr = k_ppr_0.*exp(-F.*L_int./Vt);
                %test = f_of_zeta(pi/2, F)   %f of zeta function works
                zeta_eq = Rtotal./k_ppd_eq;   %use ./ to tell it that these are not matrixes! otherwise get error.
                zeta = (Jx./L_int+ k_ppr.*zeta_eq + Rtotal)./(k_ppr+k_ppd);
                
                Un(l) = k_ppd.*zeta - Rtotal;
                Up(l-1) = k_ppd.*zeta - Rtotal;
                if(manual_Rp)
                    temp = ( E_gap - F*L_int+ Vt*log(n_full(l+1)*p_full(l)*N^2/(N_LUMO*N_HOMO)) )/(q*L_int*Rp_man);  %IN TEH FIRST TRIAL, N_FULL L+1 IS NOT DEFINED!!! IS JUST = 0!!: NEED TO FIX THIS!!
                    Un(l) = (Un(l) - temp);   %similar to Setbn, l corresponds to l+1 in the full  count...
                    Up(l-1) = (Up(l-1) - temp);  %JUST THIS ELEMENT is non-zero --> corresponds to like did the Setbp in Gummel  iterations --> l-1 corresponds to l in the full count.
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Compute the functions and Jacobian
            for i = 1:num_elements   %I DEFINED FUNCTIONS TO COMPUTE 1 ELEMENT AT A TIME SO MORE EFFICIENT IN JACOBIAN!
                Fv(i) = Poissonfnc(fullV, n(i), p(i), nT(i), pT(i),  i); %need to pass it i also, b/c using entire fullV array b/c of laplacian
                Fn(i) = Continuity_n_fnc(Jn(1,i+2), Jn(1,i+1), Un(i)); %just need to pass 2 Jn values --> NOTE: Jn's definced as i.e. Jn_1+1/2 = Jn(2) is leftmost Jn
                Fp(i) = Continuity_p_fnc(Jp(1,i+2), Jp(1,i+1), Up(i));  %NOTE:I need to pass Up(i) b/c need to define Up array--> non-zero only at 1 spot...
                
                %             FnT(i) = nTraps_fnc(n(i), nT(i));
                %             FpT(i) = pTraps_fnc(p(i), pT(i));
            end
            %Fu is upudated outside of the for loop, b/c has only 2 elements
            %that are non-zero
            %         Fu = Net_Genfnc(U, k_ppd, zeta, Rtotal);
            
            
            %NOTE: WHEN DOING THE ABOVE CALLS, MOBILITIES ARE STILL SCALED...
            %--> are mob/mobil....
            
            
            F_combined = [Fv, Fn, Fp];
            old_F_combined = F_combined;
            RHS = -F_combined;
            
            Jac = FindJacobian(V, n, p, fullV, n_full, p_full, Un, Up, B_n, B_p);  %in FindJacobian will call the functions to do the computations BUT also need the values of variables,
            % don't need nT and pT values b/c no numerical
            %derivatives computed from those.
            %NOTE: Un and Up are passed but are currently NOT used
            
            
            update_vector = Jac\RHS.'; %need to take transpose of RHS for matrices to have proper dimensions...
            
            %damped update
            old_solution = solution;
            solution = solution + w_newton*update_vector;
            
            %now update/redefine the individual solutions so can use in above
            %equations...
            V = solution(1:num_elements).';
            n = solution(num_elements+1:2*num_elements).';
            p = solution(2*num_elements+1:3*num_elements).';
            %         U = solution(3*num_elements+1:4*num_elements).';
            %         nT = solution(3*num_elements+1:4*num_elements).';
            %         pT = solution(4*num_elements+1:5*num_elements).';
            for i = 1:num_elements
                if(p(i) <0.0)
                    p(i) = 0;
                    negative_nps = true;
                end
                if(n(i) <0.0)
                    n(i) = 0;
                    negative_nps = true;  %will allow to ignore error associated with setting these to 0
                end
            end
            
            
            %for convergence analysis
            p_solution(iter) = p(20);    % save p at random point for each iter for convergence analysis
            
            fullV = [fullV(1), V, fullV(num_cell+1)];
            n_full = [n_full(1), n, n_full(num_cell+1)];
            p_full = [p_full(1), p, p_full(num_cell+1)];
            %        nT_full = [nT_full(1), nT, nT_full(num_cell+1)];
            %        pT_full = [pT_full(1), pT, pT_full(num_cell+1)];
            
            
            %Recompute F_combined...
            B_n = BernoulliFnc_n(fullV);
            B_p = BernoulliFnc_p(fullV);
            
            % Calculate drift diffusion J's
            % Use the SG definition--> BUT PULL OUT THE COEFFS (and q's
            % cancel), and coeff's go to Cn and Cp
            for i = 1:num_cell        %verified that this is right
                %NOTE: NO coefficients, mobilities here, b/c those were moved
                %over to the Cn --> these J's are passed to Continuity  fnc's:
                %which use Jr - Jl +CnUn = 0 etc..., so coeffs are in the Cn
                Jp(1,i+1) = -(p_full(i+1)*B_p(2,i+1)-p_full(i)*B_p(1,i+1));      %(1,i)  b/c want a row vector  %need an N b/c my p's are scaled by N. Extra +1's are b/c B starts at i=2
                Jn(1,i+1) =  (n_full(i+1)*B_n(1,i+1)-n_full(i)*B_n(2,i+1));
            end
            
            
            %recompute F:
            for i = 1:num_elements   %I DEFINED FUNCTIONS TO COMPUTE 1 ELEMENT AT A TIME SO MORE EFFICIENT IN JACOBIAN!
                Fv(i) = Poissonfnc(fullV, n(i), p(i), nT(i), pT(i),  i); %need to pass it i also, b/c using entire fullV array b/c of laplacian
                Fn(i) = Continuity_n_fnc(Jn(1,i+2), Jn(1,i+1), Un(i)); %just need to pass 2 Jn values --> NOTE: Jn's definced as i.e. Jn_1+1/2 = Jn(2) is leftmost Jn
                Fp(i) = Continuity_p_fnc(Jp(1,i+2), Jp(1,i+1), Up(i));  %NOTE:I need to pass Up(i) b/c need to define Up array--> non-zero only at 1 spot...
            end
            %             Fu = Net_Genfnc(U, k_ppd, zeta, Rtotal);
            
            F_combined = [Fv, Fn, Fp];
            
            
            
            old_error = error_np;
            %maybe base error on how far F_combined is from 0!!!
            
            %when we do error, we must IGNORE any elements where n and p was
            %set to 0.
            %minority n's are located at right 1/2 of device
            %minority p's are at left half of device:
            %recall F_combined is defined as [Fv, Fn, Fp]..
            %so we need to ignore the 1st 1/2 of Fn's and the 2nd 1/2 of Fp's
            error_array = F_combined;
            if(negative_nps ==true)
                clear error_array
                error_array =  [Fn(l+2:num_elements), Fp(1:l-1)];   %I BELIEVE the problem is with the Fn(50) and Fn(51)... there the errors remain high!
                %Fp's seem fine --> all errors are e-17 level!..., maybe don't
                %even need to exclude part of the array
                %Also i think Fv needs to be excluded --> errors are high--> bc
                %when set things to 0 --> is affecting the n-p term in Fv...
            end
            
            error_np = max(abs(error_array))   %update_vector should be 0 when converged
            
            %store error for convergence plotting
            error_np_array(iter) = error_np;
            iterations(iter) = iter;
            
            %if(Va == Va_max) %only for last run
            p_solution(iter) = p(20);    % save p at random point for each iter for convergence analysis
            %end
            %------------------------------------------------------------------------------------------------
            
            
            %auto relax tolerance and Newton damping if not converging
            if(error_np>= old_error)
                not_cnv_cnt = not_cnv_cnt+1;
            end
            if(not_cnv_cnt>200)
                w_newton = w_newton/2.;
                %             Newton_tolerance = Newton_tolerance*10;
                not_cnv_cnt = 0;  %reset the count
            end
            
            
            w_newton
            Newton_tolerance
            
            t =  toc
            
            %print to file (opened before Newton loop started) an element of
            %update vector...
            fprintf(fupdate,'%.8e %.8e %.8e %.8e \r\n ',Newton_iter, update_vector(40), max(abs(F_combined)), t); %NOTE: will print out all positives for F_combined...
            
            
            iter
            iter =  iter+1;
            
            
        end
    end
    
    fclose(fupdate);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Post processing once a Va is converged
    %clear the previous Jp's and Jn's
    
    clear Jp Jn
    
    % Calculate drift diffusion J's
    % Use the SG definition
    for i = 1:num_cell-1        %verified that this is right
        %HERE WE keep the mobilities, b/c these J's are for the JV curve
        %--> want the true J values here...
        Jp_temp(1,i) = -(q*p_mob*Vt*N/dx)*(p_full(i+1)*B_p(2,i+1)-p_full(i)*B_p(1,i+1));      %(1,i)  b/c want a row vector  %need an N b/c my p's are scaled by N. Extra +1's are b/c B starts at i=2
        Jn_temp(1,i) =  (q*n_mob*Vt*N/dx)*(n_full(i+1)*B_n(1,i+1)-n_full(i)*B_n(2,i+1));
    end
    
    for i =  2:num_cell
        Jp(1,i) = Jp_temp(1,i-1);     %define Jp as on the right side (i.e. i+1/2 goes to i+1).
        Jn(1,i) = Jn_temp(1,i-1);
    end
    
    %Setup for JV curve
    if(Va_cnt>0)
        V_values(Va_cnt,:) = Va;
        Jp_all (Va_cnt,:) = Jp;  %is writing each Jp in a row corresponding to Va_cnt
        Jn_all(Va_cnt,:) = Jn;
        J_total(Va_cnt,:) =  Jp + Jn;
    end
    %RIGHT NOW SOMEHTING WRONG WITH  ABOVE SETUP SO COMMENT  OUT FOR NOW....
    
    %numerical E result calculation
    for i = 2:num_cell+1
        E(i) = -Vt*(fullV(i) - fullV(i-1))/dx;  %*Vt to rescale  back to normal units: RECALL THAT I DEFINED dV as just V(i+1) - V(i) and here we need dV/dx!!
    end
    
    %Save data
    str = sprintf('%.8f',Va);
    filename = [str 'V.txt']
    fid = fopen(fullfile('C:\Users\Tim\Documents\Duxbury_group_research\Drift_Diffusion_OLED\Newtons_method_OLED\results_output',filename),'w');   %w: Open or create new file for writing
    
    if(Va_cnt ==0)
        equil = fopen(fullfile('C:\Users\Tim\Documents\Duxbury_group_research\Drift_Diffusion_OLED\Newtons_method_OLED\results_output\Equil.txt'),'w');   %w: Open or create new file for writing
        for i = 2:num_cell
            fprintf(equil,'%.8e %.8e %.8e %.8e %.8e\r\n ',x(i), Vt*fullV(i), N*p_full(i), N*n_full(i), E(i));
            %f means scientific notation, \r\n: start new line for each value
            %for each new desired column, put its own % sign
        end
        fclose(equil);
    end
    
    if(Newton_iter == 0) %then converged fully with Gummel iter and Un is just a value
        Rate = Un;
    else
        Rate = Un(l); %if converged with newton, pick out the element of rate at interface
    end
    
    if(Va_cnt > 0)
        for i = 2:num_cell
            fprintf(fid,'%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\r\n ',x(i), Vt*fullV(i), N*p_full(i), N*n_full(i), J_total(Va_cnt,i), E(i), Rate, w_newton, Newton_tolerance);
            %f means scientific notation, \r\n: start new line for each
            %value  %NEED OT USE Un(i-1) in above,(i-1) for matrix dimensions to match b/c otherwise, Un is a
            %vector and messes up the output file
            %for each new desired column, put its own % sign
            
            %NOTE THE PRINT OUT OF uN  STILL HAS ISSUES B/C IF DOESN'T GO
            %THROUGH NEWTON, THEN WILL PRINT ALL 0'S HERE!! --> NEED TO DO
            %IF STATEMENTS AND SOMETHING LIEK: IF USE NEWTON, THEN Un(i-1)
            %is used, if Gummel, then Un  is used...
            % just i.e. rate =  Un or rate =  Un(i-1)
        end
    end
    fclose(fid);
    
end
toc


%% JV setup:  trying  to print all J values to file didn't work, b/c there's
% too many J values!
file2 = fopen(fullfile('C:\Users\Tim\Documents\Duxbury_group_research\Drift_Diffusion_OLED\Newtons_method_OLED\results_output\JV.txt'),'w');
% %fullfile allows to make filename from parts

for i = 1:Va_cnt
    fprintf(file2, '%.8e %.8e \r\n', V_values(i,1), J_total(i,l-4));
end

fclose(file2);


%% Final Plots: done for the last Va

str = sprintf('%.4g', Va);

h1 = plot(x,N*p_full);
hold on
title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);
axis([-inf inf 0 inf]);

figure
plot(x,log(N*p_full));
hold on
%title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
%xlabel('Position ($m$)','interpreter','latex','FontSize',14);
%ylabel({'Log of hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);
%axis([-inf inf 0 inf]);
plot(x,log(N*n_full));
hold on
title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Log of carrier densities ($1/m^3$)'},'interpreter','latex','FontSize',14);
axis([-inf inf 0 inf]);

figure
plot(x,N*n_full);
hold on
title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Electron density ($1/m^3$)'},'interpreter','latex','FontSize',14);
axis([-inf inf 0 inf]);

%Plot E
figure
plot(x(2:num_cell), E(2:num_cell))      %We don't plot left bndry pt., b/c E there is not calculated. dV starts at i=2.
hold on
title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);
if -(fullV(2)-fullV(1)) < 0.0
    y_min = -inf;  %allow neg. y-min if necessary
else
    y_min = 0;
end
axis([-inf inf y_min inf]);


%Plot potential (fullV)
figure
plot(x, fullV)
hold on
title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Electric Potential (V)'},'interpreter','latex','FontSize',14);


%  figure;
%  h3 = plot(x(2:num_cell),Jp(2:num_cell));
%  hold on
%  title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
%  xlabel('Position ($m$)','interpreter','latex','FontSize',14);
%  ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);

% JV curve
figure
h4 = plot(V_values,J_total(:,l-4));
hold on
xlabel('Voltage (V)','interpreter','latex','FontSize',14);
ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);

%convergence analysis
figure
h5 = plot(iterations, p_solution);
title(['p convergence', str, 'V'],'interpreter','latex','FontSize',24);
xlabel('Iterations','interpreter','latex','FontSize',24);
ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',24);

figure
h6 = plot(iterations,error_np_array);
title(['Newton method convergence for', str, 'V'],'interpreter','latex','FontSize',24);
xlabel('Iterations','interpreter','latex','FontSize',24);
ylabel({'Error'},'interpreter','latex','FontSize',24);

hold off
