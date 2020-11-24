%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Solving Poisson + Drift Diffusion eqns using
%                   Scharfetter-Gummel discretization
%
%     This models a perovskite device: PEDOT:PSS/perovskite/C60/BCP
%      The HTL is PEDOT:PSS and the ETL has 2 layers (C60 and BCP)
%
%                     Timofey Golubev (2017.11.10)
%             NOTE: i=1 corresponds to x=0, i=nx to x=L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;   

%% Declare global variables  (accessible by all functions)
global L num_cell num_elements Cp Cn  N q kb T Vt dx p_mob n_mob l_HTL_int;
global l_ETL_int HTL_traps ETL_traps l_BCP_int BCP_int_VBstep BCP_int_CBstep;
global N_LUMO N_HOMO cap_n cap_p k_rec phi_c phi_a n1 p1 HTL_int_VBstep ETL_int_VBstep ETL_int_CBstep HTL_int_CBstep G_max;

%% Parameters
%Thicknesses (in m)
L_HTL = 1.5*10^-9;  %PEDOT:PSS 
L_perovskite = 320*10^-9;
L_ETL = 20.0*10^-9;  %C60
L_BCP = 7.5*10^-9;
L = (L_HTL+L_perovskite+L_ETL + L_BCP);   %device length in meters: 
dx = 0.5*10^-9;

num_cell = L/dx;            % number of cells   
l_HTL_int = L_HTL/dx + 1;                           % position of HTL/perovskite interface (+1 b/c position 1 corresponds to x = 0)   
l_ETL_int = (L_HTL + L_perovskite)/dx + 1;          % position of perovskite/C60 (ETL) interface 
l_BCP_int = (L_HTL + L_perovskite + L_ETL)/dx + 1;  % position of C60/BCP (ETL) interface

% Physical Constants
q =  1.60217646*10^-19;         % elementary charge, C
kb = 1.3806503*10^-23;          % Boltzmann const., J/k
T = 296.;                       % temperature (K)
epsilon_0 =  8.85418782*10^-12; % F/m
Vt = (kb*T)/q;       % thermal voltage
N_LUMO = 8.1*10^24;  % density of states
N_HOMO = 8.1*10^24;
N = 8.1*10^24.;      % scaling factor helps CV be on order of 1 for improved numerical stability
  
%Energy levels (in eV)
perov_CB = -3.9;
perov_VB = -5.4;   
ETL_LUMO = -4.2;
ETL_HOMO = -5.9;  
HTL_HOMO = -5.1;  
HTL_LUMO = -3.5;
BCP_LUMO = -4.32;
BCP_HOMO = -7.0;
E_gap = 1.5;        % Perovskite bandgap

phi_a = 0.3;	% Anode injection barrier
phi_c = 0.0;	% Cathode injection barrier
G_max = 1.97*10^28;  % scaling for generation rate

%Position dependent dielectric constants
epsilon(1:l_HTL_int) = 3.0*epsilon_0;                % HTL
epsilon(l_HTL_int +1: l_ETL_int) = 24.1*epsilon_0;  % perovskite  
epsilon(l_ETL_int+1:l_BCP_int) = 4.25*epsilon_0;    % ETL (C60)
epsilon(l_BCP_int+1:num_cell+1) = 4.25*epsilon_0;   % ETL (BCP)

% Holes
p_mob(1:l_HTL_int) = 4.5*10^-6;           % HTL          
p_mob(l_HTL_int + 1: l_ETL_int) = 10^-4;  % perovskite  
p_mob(l_ETL_int+1:l_BCP_int) = 1.6*10^-4; % ETL (C60)
p_mob(l_BCP_int+1:num_cell +1)= 2.*10^-9; % ETL (BCP)
% Electrons
n_mob(1:l_HTL_int) = 4.5*10^-6;           % HTL
n_mob(l_HTL_int + 1: l_ETL_int) = 10^-4;  % perovskite  
n_mob(l_ETL_int+1:l_BCP_int) = 1.6*10^-4; % ETL (C60)
n_mob(l_BCP_int+1:num_cell +1)= 2.*10^-9; % ETL (BCP)
mobil =  5.0*10^-6;       %scaling for mobility (for numerical stability)

Vbi = 4.8-4.3; % built-in voltage

% calculated energy barriers
HTL_int_VBstep =  abs(perov_VB - HTL_HOMO);     
HTL_int_CBstep = abs(perov_CB - HTL_LUMO);
ETL_int_VBstep = abs(perov_VB - ETL_HOMO);
ETL_int_CBstep = abs(perov_CB - ETL_LUMO);
BCP_int_VBstep = abs(perov_VB - BCP_HOMO);      
BCP_int_CBstep = (ETL_LUMO - BCP_LUMO); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Voltage step setup (JV curve start-stop)
Va_min = -0.01;              % start of JV curve
Va_max = 1.3;                % stop of JV curve
increment = 0.01;            % voltage step   
num_V = floor((Va_max-Va_min)/increment)+1; % number of steps

%Simulation parameters
w_i = 0.2;                   % set up of weighting factor (relaxation for stability)
tolerance = 5*10^-12;        % error tolerance  
tolerance_i = 5*10^-12;      % initial error tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Recombination parameters: 
k_rec = 6*10^-17; % Langevin recombination coefficient

%SRH trap-assisted recombination 
cap_n = 10^-13;  
cap_p = 10^-13;

%Shockley-Read-Hall (SRH) recombo assumed to occur in 2 mesh cells at interface (so 1nm region)
HTL_traps = 5*10^21; 
ETL_traps = 6*10^20;

E_trap = perov_VB + E_gap/2.0;  % Assume that energy of the traps are in middle of bandgap (this is where recombination is most effective)
n1 = N_LUMO*exp(-(perov_CB - E_trap)/Vt);  
p1 = N_HOMO*exp(-(E_trap - perov_VB)/Vt);

%% Domain Discretization
a_1=0; a_2=L; x=linspace(a_1,a_2,num_cell+1);  %x is positions array, +1 necessary b/c linspace needs (# of pts = # of cells +1)  
nx = length(x);  
num_elements = nx-2;

%% Boundary Conditions
n_full(1) = N_LUMO*exp(-(E_gap - phi_a)/Vt);       % anode 
p_full(1) = N_HOMO*exp(-phi_a/Vt);                
n_full(num_cell+1) = N_LUMO*exp(-phi_c/Vt);        % cathode
p_full(num_cell+1) = N_HOMO*exp(-(E_gap - phi_c)/Vt);

% scaling
n_full = n_full/N;
p_full = p_full/N;

%% Initial Conditions
min_dense = min(n_full(1),  p_full(num_cell+1)); 

% load generation rate input
gen_rate_input =load('gen_rate.txt');

for i = 1:num_elements         
    p(i) = min_dense; 
    n(i) = min_dense;  
end

%set up initial guess of V make V linearly decreasing     
fullV(1)          = -Vbi/(2*Vt);  
fullV(num_cell+1) = Vbi/(2*Vt); 
diff = (fullV(num_cell+1)-fullV(1))/num_cell;   
for i = 2:num_elements+1
    fullV(i) = fullV(i-1)+diff;
end
V_initial = fullV(2:num_elements+1);  %shift to correspond to V in matrix: this will be used as initial guess in QR solver
V_initial = V_initial.'; 


%% Main voltage loop
Va_cnt = 0;
for Va_cnt = 0:num_V+1   %+1 b/c 1st Va is the equil. Va = 0 run
    not_converged = false;
    not_cnv_cnt = 0;
    
    %stop the calculation if tolerance becomes too high
    if(tolerance >10^-4)
        break
    end
     
    %allocate
    V = zeros(nx-2,1);  
    
    if(Va_cnt ==0) 
        % Convergence from an initial guess on the first voltage step is
        % difficult. Therefore, we relax the tolerance, and use a very 
        % small relaxation (w).
        tolerance = tolerance*10^4;       
        w = 0.01;   
        Va = 0;
        Un = zeros(num_elements,1);
        Up = zeros(num_elements,1);
        V = V_initial;
    end
    if(Va_cnt ==1)
        tolerance = tolerance_i; %reset tolerance back to initial tolerance
        w = w_i;  
    end

    if(Va_cnt >0)
        Va = Va_min+increment*(Va_cnt-1);     %set Va value
    end
     % apply voltage boundary conditions
     fullV(1)          = -(Vbi-Va)/(2*Vt);  
     fullV(num_cell+1) = (Vbi-Va)/(2*Vt);
     
     % Concatenate charge densities BCs with charge densities solution
     p_full = [p_full(1), p, p_full(num_cell+1)]; 
     n_full = [n_full(1), n, n_full(num_cell+1)];

    %timing
    tic
    
    %% Setup the solver
    iter = 1;
    error_np =  1.0;
    
    % set up AV matrix: these never change so do outside of loop
    
    % we use epsilon/epsilon_0 b/c want to have only the relative
    % permitivities as coeffs in matrix so that it is better-scaled
    AV_val = SetAV_val(epsilon/epsilon_0);  
    bV = zeros(num_elements,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Solver Loop         
    Cp = dx^2/(Vt*N*mobil);          % scaled p_mob and n_mob are inside matrices
    Cn = dx^2/(Vt*N*mobil);  
    CV = N*dx^2*q/(epsilon_0*Vt);    % relative permitivity was moved into the matrix

    while error_np > tolerance

        %% Poisson equation with tridiagonal solver
        % setup bV
        for i = 1: num_elements  
           bV(i,1) = CV*((n(i))-(p(i)));   
        end

        %for bndrys 
        bV(1,1) = bV(1,1) - (epsilon(1)/epsilon_0)*fullV(1);        
        bV(num_elements,1) = bV(num_elements,1) - (epsilon(num_cell+1)/epsilon_0)*fullV(num_cell+1);
        
        oldV = V;
        newV = AV_val\bV;

        if(iter >1)
          V = newV*w + oldV*(1.-w);
        else
          V = newV;  % no mixing for 1st iter
        end

        % Concatenate V solution with bndry pts
        fullV = [fullV(1); V; fullV(num_cell+1)];  % note: bndry values are already defined above
        fullV = fullV.';              
%------------------------------------------------------------------------------------------------        
        %% Calculate recombination rates  
        if(Va_cnt > 0)
             G = GenerationRate(gen_rate_input);
        else
            G = zeros(num_elements,1);
        end

         R_Langevin_array = R_Langevin(n_full,p_full); 
         R_SRH_HTL_array = R_SRH_HTL(n_full, p_full);
         R_SRH_ETL_array = R_SRH_ETL(n_full, p_full);

         % Net carrier generation equals generation rate minus
         % recombination rate.
         for i = 1:num_elements    
            Un(i) = G(i)- R_Langevin_array(i);
         end
         Up = Un;  % assume holes generation rate equals that of electrons
         
         % 5 nm trap-assited recombination regions, near the ETL and HTL
         % boundaries.
         for i =  l_ETL_int-9:l_ETL_int  
             Up(i-1) = Up(i-1) - R_SRH_ETL_array(i);       
             Un(i-1) = Un(i-1)- R_SRH_ETL_array(i);
         end
       
         for i = l_HTL_int+1:l_HTL_int+10      
              Un(i-1) = Un(i-1) - R_SRH_HTL_array(i);
              Up(i-1) = Up(i-1)- R_SRH_HTL_array(i);
         end
       
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Now solve eqn for n and p  
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
        newp = p_sol.'; 
        
        old_n = n;
        n_sol = An_val\bn;
        newn = n_sol.';
        
        % If get negative p's or n's, make them equal 0
        % the negatives can arise during iterating due to numerical 
        % oscillations
        for i = 1:num_elements        
            if(newp(i) <0.0)
                newp(i) = 0;
            end    
            if(newn(i) <0.0)
                newn(i) = 0;
            end
        end
        
        % calculation of error--> iteration-to-iteration change
        % in charge densities. Once error is less than the tolerance,
        % we have found a converged solution.
        old_error =  error_np;
        count = 0;
        error_np_matrix = zeros(1,num_elements); 
        for i = 1:num_elements
            if(newp(i) ~=0 && newn(i) ~=0)
                count = count+1;  %counts number of non zero error calculations
                error_np_matrix(count) = (abs(newp(i)-old_p(i))+abs(newn(i)-old_n(i)))/abs(old_p(i)+old_n(i)); 
            end
        end
        error_np = max(error_np_matrix) 
          
        % auto decrease w (relaxation factor) if solution is not converging
        if(error_np>= old_error)
            not_cnv_cnt = not_cnv_cnt+1;
        end
        if(not_cnv_cnt>1000)
            w = w/2.;
            tolerance = tolerance*10;
            not_cnv_cnt = 0;  %reset the count
        end
        w
        tolerance
        
        %weighting (relaxation) for numerical stability
        p = newp*w + old_p*(1.-w);
        n = newn*w + old_n*(1.-w);
        
        iter =  iter+1    
          
       %unscale mobilities
       p_mob = p_mob*mobil;
       n_mob = n_mob*mobil;
    
       p_full = [p_full(1), p, p_full(num_cell+1)]; % bndrys already defined above
       n_full = [n_full(1), n, n_full(num_cell+1)];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % Calculate drift diffusion J's
    % Use the SG definition
    for i = 1:num_cell-1        
        Jp_temp(1,i) = -(q*Vt*N/dx)*(p_mob(i+1)*p_full(i+1)*B_p(2,i+1)-p_mob(i+1)*p_full(i)*B_p(1,i+1));     
        Jn_temp(1,i) =  (q*Vt*N/dx)*(n_mob(i+1)*n_full(i+1)*B_n(1,i+1)-n_mob(i+1)*n_full(i)*B_n(2,i+1)); 
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
        E(i) = -Vt*(fullV(i) - fullV(i-1))/dx;  %*Vt to rescale  back to normal units: reccall that dV is just V(i+1) - V(i) and here we need dV/dx
    end
    
    % Save data
    str = sprintf('%.2f',Va);
    filename = [str 'V.txt'] 
    fid = fopen(fullfile(filename),'w');   
  
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
toc
 
%% JV setup:  
 file2 = fopen(fullfile('JV.txt'),'w');

for i = 1:Va_cnt
   fprintf(file2, '%.8e %.8e \r\n', V_values(i,1), J_total(i,floor(num_cell/2)));
end

fclose(file2);

%% Final Plots: done for the last Va

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
 plot(x(2:num_cell), E(2:num_cell))      % We don't plot left bndry pt., b/c E there is not calculated. dV starts at i=2.
 hold on
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);
 if -(fullV(2)-fullV(1)) < 0.0
     y_min = -inf;  % allow neg. y-min if necessary
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

% JV curve
 figure
 h4 = plot(V_values,J_total(:,floor(num_cell/2)));  %take J from middle of pervoskite
 hold on
 xlabel('Voltage (V)','interpreter','latex','FontSize',14);
 ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
 hold off
