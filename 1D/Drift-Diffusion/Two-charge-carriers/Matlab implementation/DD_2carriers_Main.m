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
clear; close all; clc;     %clear improves performance over clear all, and still clears all variables.

global G_max num_cell num_elements n_mob p_mob Cp Cn Vt

%% Parameters
L = 10e-9;      %device length in meters
dx = 1e-9;                        %mesh size
num_cell = floor(L/dx);

%Electronic density of states of holes and electrons
N_VB = 10^24;         %density of states in valence band (holes)
N_CB = 10^24;         %density of states in conduction bands (electrons)
E_gap = 1.5;          %bandgap (in eV)

N = 10^24.;            %scaling factor helps CV be on order of 1

% Physical Constants
q =  1.60217646*10^-19;         %elementary charge, C
kb = 1.3806503*10^-23;          %Boltzmann const., J/k
T = 296.;                       %temperature
epsilon_0 = 8.85418782*10^-12;     %F/m
Vt = (kb*T)/q;

%injection barriers
inj_a = 0.2;	%at anode
inj_c = 0.1;	%at cathode

%work functions of anode and cathode
WF_anode = 4.8;     
WF_cathode = 3.7;

Vbi = WF_anode - WF_cathode +inj_a +inj_c;  %built-in field

G_max = 4*10^27; %for photogeneration rate normalization

% Dielectric constants (can be made position dependent, as long as are
% piecewise-constant)
epsilon(1:num_cell+1) = 3.0*epsilon_0; 

% Carrier mobilities (can be made position dependent, as long as are
% piecewise-constant)
p_mob(1:num_cell+1) = 4.5*10^-6;       
n_mob(1:num_cell +1) = 4.5*10^-6;      

mobil = 5.*10^-6;                    %scaling for mobility (for numerical stability when solving matrix equations)

%% JV sweep parameters
Va_min = -0.5;               %volts
Va_max = -0.45;
increment = 0.01;         %for increasing V
num_V = floor((Va_max-Va_min)/increment)+1;   %number of V points

%Simulation parameters
w_eq = 0.01;               %For the 1st iteration (equilibrium run) to converge need a small weighting factor
w_i = 0.2;                 %set up of weighting factor
tolerance = 5*10^-12;        %error tolerance
tolerance_i =  5*10^-12;

%% Domain Discretization
a=0; b=L; x=linspace(a,b,num_cell+1); dx=(b-a)/num_cell;   %x is positions array, +1 necessary b/c linspace needs (# of pts = # of cells +1)
nx = length(x);
num_elements = nx-2;

%% Boundary Conditions
n_full(1) = N_CB*exp(-(E_gap-inj_a)/Vt);
p_full(1) = N_VB*exp(-inj_a/Vt);
n_full(num_cell+1) = N_CB*exp(-inj_c/Vt);
p_full(num_cell+1) = N_VB*exp(-(E_gap-inj_c)/Vt);

n_full = n_full/N;
p_full = p_full/N;

%% Initial Conditions
min_dense = min(n_full(1),  p_full(num_cell+1));
for i = 1:num_elements
    p(i) = min_dense;
    n(i) = min_dense;
end

%set up initial guess for electric potential V (simplest is just linear)
fullV(1) = -((Vbi)/(2*Vt)-inj_a/Vt);
fullV(num_cell+1) = (Vbi)/(2*Vt)-inj_c/Vt;
diff = (fullV(num_cell+1)-fullV(1))/num_cell;
for i = 2:num_elements+1
    fullV(i) = fullV(i-1)+diff;
end
V_initial = fullV(2:num_elements+1);  %shift to correspond to V in  matrix: this will be used as initial guess in QR solver
V_initial = V_initial.';   %transpose to column vector

V = zeros(nx-2,1);  %allocate array for Poisson solution

%% Main voltage loop
Va_cnt = 0;
for Va_cnt = 0:num_V +1
    not_converged = false;
    not_cnv_cnt = 0;
    
    %stop the calculation if tolerance becomes too high
    if(tolerance >10^-5)
        break
    end
    
    %1st iteration is to find the equilibrium values (no generation rate)
    if(Va_cnt ==0)
        tolerance = tolerance*10^2;       %relax tolerance for equil convergence
        w = w_eq;   %set w to be small for equilibrium convergence
        Va = 0;
        Un = zeros(num_elements,1);
        Up = zeros(num_elements,1);
        V = V_initial;
    end
    
    if(Va_cnt ==1)
        tolerance = tolerance_i;       %reset tolerance back
        w=w_i;
    end
    
    if(Va_cnt >0)
        Va = Va_min+increment*(Va_cnt-1);     %set Va value
        %Va = Va_max-increment*(Va_cnt-1);    %decrease Va by increment in each iteration
    end
    
    %Boundary conditions
    fullV(1) = -((Vbi  -Va)/(2*Vt)-inj_a/Vt);
    fullV(num_cell+1) = (Vbi- Va)/(2*Vt) - inj_c/Vt;
    
    %% Setup the solver
    iter = 1;
    error_np =  1.0;
    
    %set up AV: these never change so do outside of while loop
    AV_val = SetAV_val(epsilon/epsilon_0);  %NOTE: we use epsilon/epsilon_0 b/c want to have only the relative permitivities as coeffis in matrix
    bV = zeros(num_elements,1);
    
    %% Solver Loop
    %Scaling coefficients
    Cp = dx^2/(Vt*N*mobil);          %note: scaled p_mob and n_mob are inside matrices
    Cn = dx^2/(Vt*N*mobil);
    CV = N*dx^2*q/(epsilon_0*Vt);    %relative permitivity was moved into the matrix
    
    while error_np > tolerance
        
        %Poisson equation with tridiagonal solver
        % setup RHS
        for i = 1: num_elements
            bV(i,1) = CV*(n(i)-p(i));
        end
        %for bndrys
        bV(1,1) = bV(1,1) - (epsilon(1)/epsilon_0)*fullV(1);
        bV(num_elements,1) = bV(num_elements,1) - (epsilon(num_cell+1)/epsilon_0)*fullV(num_cell+1);
        
        oldV = V;
        newV = AV_val\bV;
        
        %Mix V's: I think I must to prevent mixing on 1st iter b/c on 1st iter,
        %need to ensure that V changes according to Va, not get mixed with the old Va
        if(iter >1)
            V = newV*w + oldV*(1.-w);
        else
            V = newV;  %no mixing for 1st iter
        end
        
        %make V with bndry pts
        fullV = [fullV(1); V; fullV(num_cell+1)];       %note: bndry values are already defined above
        fullV = fullV.';             %transpose
        
        % Calculate net carrier generation rates
        if(Va_cnt > 0)
            G = GenerationRate();  
        else
            G = zeros(num_elements,1);
        end
        
        %recombination
        %here can call a function to calculate recombination rates. In this
        %example I just set it to 0
        Recombination(1:num_elements) = 0;
        
        for i = 1:num_elements    %NOTE: i =1 corresponds to x = 2nm, and Langevin array i's are based on real x values, so need to use i+1 in langevin array
            Un(i) = G(i)- Recombination(i);
        end
        Up = Un;
        %------------------------------------------------------------------------------------------------
        %% now solve eqns for electron and hole densities (n and p respectively)
        %scale mobilities
        p_mob = p_mob/mobil;
        n_mob = n_mob/mobil;
        
        Bp = BernoulliFnc_p(fullV);
        Bn = BernoulliFnc_n(fullV);
        Ap_val = SetAp_val(Bp);
        An_val = SetAn_val(Bn);
        bp = Setbp(Bp, p_full, Up);
        bn = Setbn(Bn, n_full, Un);
        
        old_p = p;
        p_sol = Ap_val\bp;
        newp = p_sol.';    %transpose
        
        old_n = n;
        n_sol = An_val\bn;
        newn = n_sol.';
        
        % if get negative p's or n's, make them equal 0
        for i = 1:num_elements
            if(newp(i) <0.0)
                newp(i) = 0;
            end
            if(newn(i) <0.0)
                newn(i) = 0;
            end
        end
        
        old_error = error_np;
        count = 0;
        error_np_matrix = zeros(1,num_elements); %need to reset the matrix b/c when more newp's become 0, need to remove that error element from matrix.
        for i = 1:num_elements
            if(newp(i) ~=0 && newn(i) ~=0)
                count = count+1;  %counts number of non zero error calculations
                error_np_matrix(count) = (abs(newp(i)-old_p(i))+abs(newn(i)-old_n(i)))/abs(old_p(i)+old_n(i));
            end
        end
        error_np = max(error_np_matrix)
        
        %auto decrease w if not converging
        if(error_np>= old_error)
            not_cnv_cnt = not_cnv_cnt+1;
        end
        if(not_cnv_cnt>2000)
            w = w/2.;
            tolerance = tolerance*10;
            not_cnv_cnt = 0;  %reset the count
        end
        
        %output w and tolerance
        w
        tolerance
        
        %weighting
        p = newp*w + old_p*(1.-w);
        n = newn*w + old_n*(1.-w);
        
        %unscale mobilities
        p_mob = p_mob*mobil;
        n_mob = n_mob*mobil;
        
        p_full = [p_full(1), p, p_full(num_cell+1)]; %bndrys already defined above
        n_full = [n_full(1), n, n_full(num_cell+1)];
     
      iter =  iter+1
       
    end
   
   
    % Calculate drift diffusion currents
    % Use the SG definition
    for i = 1:num_cell-1
         Jp_temp(1,i) = -(q*Vt*N/dx)*(p_mob(i+1)*p_full(i+1)*Bp(2,i+1)-p_mob(i+1)*p_full(i)*Bp(1,i+1));
         Jn_temp(1,i) =  (q*Vt*N/dx)*(n_mob(i+1)*n_full(i+1)*Bn(1,i+1)-n_mob(i+1)*n_full(i)*Bn(2,i+1));
% Jp_temp(1,i) = -(q*Vt*N/dx)*(p_mob(i+1)*p_full(i+1)-p_mob(i+1)*p_full(i));
%        Jn_temp(1,i) =  (q*Vt*N/dx)*(n_mob(i+1)*n_full(i+1)-n_mob(i+1)*n_full(i));
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
        E(i) = -Vt*(fullV(i) - fullV(i-1))/dx;  %*Vt to rescale  back to normal units: RECALL THAT I DEFINED dV as just V(i+1) - V(i) and here we need dV/dx
    end
    
    %Save data
    str = sprintf('%.2f',Va);
    filename = [str 'V.txt']
    fid = fopen(fullfile(filename),'w');
    %fullfile allows to make filename from parts
    
    if(Va_cnt ==0)
        equil = fopen(fullfile('Equil.txt'),'w');
        for i = 2:num_cell
            fprintf(equil,'%.8e %.8e %.8e %.8e %.8e\r\n ',x(i), Vt*fullV(i), N*p_full(i), N*n_full(i), E(i));
        end
        fclose(equil);
    end
    
    if(Va_cnt > 0)
        for i = 2:num_cell
            fprintf(fid,'%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\r\n ',x(i), Vt*fullV(i), N*p_full(i), N*n_full(i), J_total(Va_cnt,i), E(i), Un(i-1), w, tolerance);
        end
    end
    fclose(fid);
    
    
end

%% JV setup:
file2 = fopen(fullfile('JV.txt'),'w');
for i = 1:Va_cnt
    fprintf(file2, '%.8e %.8e \r\n', V_values(i,1), J_total(i,floor(num_cell/2)));
end
fclose(file2);

%% Final Plots: done for the last Va

str = sprintf('%.2g', Va);

figure
plot(x,log(N*p_full));
hold on
plot(x,log(N*n_full));
hold on
title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Log of carrier densities ($1/m^3$)'},'interpreter','latex','FontSize',14);
axis([-inf inf 0 inf]);
legend({'holes', 'electrons'});

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
plot(x, Vt*fullV)
hold on
title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Electric Potential (V)'},'interpreter','latex','FontSize',14);

%Plot Jn and Jp
figure
plot(x(2:num_elements),Jn(2:num_elements))
hold on
plot(x(2:num_elements),Jp(2:num_elements))
title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
xlabel('Position(m)','interpreter','latex','FontSize',14);
ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
legend({'Jn', 'Jp'});

% JV curve
figure
h4 = plot(V_values,J_total(:,floor(num_cell/2)));  %take J from middle of pervoskite
hold on
xlabel('Voltage (V)','interpreter','latex','FontSize',14);
ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
axis([-inf inf -inf 50])

hold off
