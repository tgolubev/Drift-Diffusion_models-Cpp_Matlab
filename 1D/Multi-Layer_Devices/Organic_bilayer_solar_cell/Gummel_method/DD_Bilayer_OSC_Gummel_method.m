%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving Poisson + Drift Diffusion eqns (for holes) using
%                   Scharfetter-Gummel discretization
%
%       This includes settings any negative newn and newp solutions to 0
%
%                  Coded by Timofey Golubev (2017.10.07)
%             NOTE: i=1 corresponds to x=0, i=nx to x=L

%RECALL THAT fullV = V/Vt: is scaled!. Same with p = p/N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;   %NOTE: clear improves performance over clear all, and still clears all variables.

%% Declare global variables  (accessible by all functions)
global L l num_cell num_elements Cp Cn  N q kb T Vt dx L_int epsilon p_mob n_mob;
global H_A H_D T_tD T_tA N_LUMO N_HOMO gap_A gap_D cap_n cap_p Eb rc k_rec phi_c phi_a gap_elec gap_hole;

%diary  %supposed to output all command window stuff to file

%% Parameters
L = 73*10^-9;             %device length in meters
l = 51;                   %position of interface
num_cell = 73;            % number of cells

p_mob = 2.0*10^-8;        %hole mobility
n_mob =  2.0*10^-7;
mobil =  2.0*10^-7;       %scaling for mobility

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

Vbi = E_gap - phi_a - phi_c;     %Built-in voltage: this controlls when current in solar cell will switch direction


N = 10^27.;                 %scaling factor helps CV be on order of 1

Va_min = 0.5;               %volts
Va_max = 2.;
increment = 0.01;            %for increasing V
num_V = floor((Va_max-Va_min)/increment)+1;

%Simulation parameters
w_i = 0.6;           %set up of weighting factor: At low V, w = 0.4 works up to about 0.52V, BUT note: this results in tolerance being reduced by 10
%and higher w for rest.
tolerance = 5*10^-13;        %error tolerance
tolerance_i = 5*10^-13;
equil_tolerance = 5*10^-12;
manual_Rp = true;

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

k_rec = q*(n_mob+ p_mob)/epsilon;      %from Langevin		%this is still correct to use n and p mob here
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

%set up initial guess of V
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
%clear iterations;
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
    
    if(Va_cnt >0)
        Va = Va_min+increment*(Va_cnt-1);     %set Va value
        %Va = Va_max-increment*(Va_cnt-1);    %decrease Va by increment in each iteration
    end
    
    %     if(Va == 0.93)
    %         diary  %diary turns diary on or off, so must only write 'diary" once !! to have it on for reminader of calculation
    %     end
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solver Loop
    Cp = dx^2/(Vt*N*mobil);   %note:now use mobil: the scaling here: b/c the actual p and n mob are inside the matrix
    Cn = dx^2/(Vt*N*mobil);
    CV = N*dx^2*q/(epsilon*Vt);
    
    %file for outputing update_vector element for convergence behavior testing
    str = sprintf('%.3f',Va);
    filename = [str 'update_vector.txt']
    fupdate = fopen(fullfile('C:\Users\Tim\Documents\Duxbury_group_research\DD-BI_Matlab_Scharfetter-Gummel\Gummel_results_output',filename),'w');   %w: Open or create new file for writing
    
    while error_np > tolerance
        
        %for above Va_cnt of 1, like in Newton, call traps only in 1st iteration, to try to
        %stabilize algo more
        if(Va_cnt > 1)
            if(iter == 1)
                [pT_full, nT_full] = CalTrap(p_full, n_full);   %calculate trap densities
                nT = nT_full(2:num_cell);   %redefine as only the points inside device
                pT = pT_full(2:num_cell);
            end
        else  %otherwise call in every iter
            [pT_full, nT_full] = CalTrap(p_full, n_full);   %calculate trap densities
            nT = nT_full(2:num_cell);   %redefine as only the points inside device
            pT = pT_full(2:num_cell);
        end
        
        %% Poisson equation with tridiagonal solver
        % setup bV
        for i = 1: num_elements
            bV(i,1) = CV*((n(i)+nT(i))-(p(i)+pT(i)));
        end
        %for bndrys
        bV(1,1) = bV(1,1) - fullV(1);
        bV(num_elements,1) = bV(num_elements,1) - fullV(num_cell+1);
        
        %call solver, solve for V
        %spparms('spumoni',2)     %allows to see details of solver
        
        oldV = V;
        newV = lsqr(AV_val,bV,10^-15,1000,[] ,[] ,V);  %can't converge to below 10^-14 on the QR here
        
        normV = norm(AV_val*V-bV)  %check convergence again
        
        %NOTE: I MUST to prevent mixing on 1st iter b/c on 1st iter,
        %need to ensure that V changes according to Va, not get mixed with the old Va
        %I checked: if mix on 1st iter also, get WRONG data --. doesn't
        %match with original version.        
        if(iter >1)
            V = newV*w + oldV*(1.-w);
        else
            V = newV;  %no mixing for 1st iter
        end
        
        %make V with bndry pts
        fullV = [fullV(1); V; fullV(num_cell+1)];       %note: bndry values are already defined above
        fullV = fullV.';   
        
        %------------------------------------------------------------------------------------------------
        %% Calculate recombination rates
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
        % below 2 lines are for \ solve
        p_sol = Ap_val\bp;
        newp = p_sol.';  %need to transpose if don't use qr
        % if using lsqr, comment above 2 lines and uncomment below
        % p_sol = lsqr(Ap_val,bp,10^-14,1000,[] ,[] ,p);
        %newp = p_sol;    %don't need to transpose b/c lsqr already gives column vector
        normp = norm(Ap_val*p_sol-bp);
        
        old_n = n;
        % below 2 lines are for \ solve
        n_sol = An_val\bn;
        newn = n_sol.';
        % if using lsqr, comment above 2 lines and uncomment below
        % n_sol = lsqr(An_val,bp,10^-14,1000,[] ,[] ,n);
        %newn = n_sol;
        normn = norm(An_val*n_sol-bn);
        
        % if get negative p's or n's, make them equal 0
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
        
        %store error for convergence plotting
        %error_np_array(iter) = error_np;
        %iterations(iter) = iter;
        
        if(Va == Va_max) %only for last run
            p_solution(iter) = p(42);    % save p at random point close to interface for each iter for convergence analysis
        end
        
        %auto decrease w if not converging
        if(error_np>= old_error)
            not_cnv_cnt = not_cnv_cnt+1;
        end
        if(not_cnv_cnt>1000)
            w = w/2.;
            %tolerance = tolerance*10;
            not_cnv_cnt = 0;  %reset the count
        end
        w
        tolerance
        
        %weighting
        p = newp*w + old_p*(1.-w);
        n = newn*w + old_n*(1.-w);
        
        iter =  iter+1;
        
        t = toc;
        
        %print to file (opened before loop started) an element of
        %update vector...
        fprintf(fupdate,'%.8e %.8e %.8e %.8e \r\n ',iter, w, tolerance, t);
        
        %unscale mobilities
        p_mob = p_mob*mobil;
        n_mob = n_mob*mobil;
        
        p_full = [p_full(1), p, p_full(num_cell+1)]; %bndrys already defined above
        n_full = [n_full(1), n, n_full(num_cell+1)];
        
    end
    fclose(fupdate);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate equil kppd so have it for future Va's
    if(Va_cnt == 0)
        F = (fullV(l+1)-fullV(l)).*Vt./L_int;
        k_ppd_eq = kppd(F);   %find equil. kppd value
    end
    
    % Calculate drift diffusion J's
    % Use the SG definition
    for i = 1:num_cell-1        %verified that this is right
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
    
    %numerical E result calculation
    for i = 2:num_cell+1
        E(i) = -Vt*(fullV(i) - fullV(i-1))/dx;  %*Vt to rescale  back to normal units: RECALL THAT I DEFINED dV as just V(i+1) - V(i) and here we need dV/dx!!
    end
    
    %Save data
    str = sprintf('%.2f',Va);
    filename = [str 'V.txt']
    fid = fopen(fullfile('C:\Users\Tim\Documents\Duxbury_group_research\DD-BI_Matlab_Scharfetter-Gummel\Gummel_results_output',filename),'w');   %w: Open or create new file for writing
    
    if(Va_cnt ==0)
        equil = fopen(fullfile('C:\Users\Tim\Documents\Duxbury_group_research\DD-BI_Matlab_Scharfetter-Gummel\Gummel_results_output\Equil.txt'),'w');   %w: Open or create new file for writing
        for i = 2:num_cell
            fprintf(equil,'%.8e %.8e %.8e %.8e %.8e\r\n ',x(i), Vt*fullV(i), N*p_full(i), N*n_full(i), E(i));
        end
        fclose(equil);
    end
    
    if(Va_cnt > 0)
        for i = 2:num_cell
            fprintf(fid,'%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\r\n ',x(i), Vt*fullV(i), N*p_full(i), N*n_full(i), J_total(Va_cnt,i), E(i), Rtotal, Un, zeta, k_ppd, k_ppr, w, tolerance);
        end
    end
    fclose(fid);
    
end
toc

%% JV setup:  trying  to print all J values to file didn't work, b/c there's
% too many J values!
file2 = fopen(fullfile('C:\Users\Tim\Documents\Duxbury_group_research\DD-BI_Matlab_Scharfetter-Gummel\Gummel_results_output\JV.txt'),'w');
for i = 1:Va_cnt
    fprintf(file2, '%.8e %.8e \r\n', V_values(i,1), J_total(i,l-4));
end
fclose(file2);

%% Final Plots: done for the last Va

matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','manual')  %THIS ENSURES THAT WHEN FORMATTING GRAPHS, it doesn't try to auto change the font sizes etc...

str = sprintf('%.2g', Va);

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

%plot on a real log plot the densities...

figure
h = semilogy(x*10^9,N*p_full,'LineWidth',1.5);
hold on
semilogy(x*10^9,N*n_full, 'LineWidth',1.5);
hold on
set(gca, 'FontSize', 32)
title(['Va =', str, 'V'],'interpreter','latex','FontSize',36);
xlabel('Position ($nm$)','interpreter','latex','FontSize',36);
ylabel({'Carrier densities ($1/m^3$)'},'interpreter','latex','FontSize',36);
xticks([0 10 20 30 40 50 60 70])
yticks([0  10 10^5 10^10 10^15 10^20 10^25])  %control placement of tick marks
%  axis([-inf inf 0 inf]);


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
%NOTE: these convergence plots are smooth unless use a Va which is close
%to where w needs to be dropped.
figure
h5 = plot(1:iter-1, p_solution);
title(['p convergence', str, 'V'],'interpreter','latex','FontSize',24);
xlabel('Iterations','interpreter','latex','FontSize',24);
ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',24);

% figure
% h6 = plot(iterations,error_np_array);
% title(['Newton method convergence for', str, 'V'],'interpreter','latex','FontSize',24);
% xlabel('Iterations','interpreter','latex','FontSize',24);
% ylabel({'Error'},'interpreter','latex','FontSize',24);

hold off
