%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving Poisson + Drift Diffusion eqns (for holes) with finite 
%           differences. Used to compare with WENO
%
%NOTE: MUCH  SMALLER MESH SIZES (THAN WENO) ARE NEEDED FOR FINITE DIFF. 
%                            METHOD TO WORK
%
%    Coded by Timofey Golubev (2017.08.11) based on original 1D wave eqn
%              code by Manuel Diaz, manuel.ade'at'gmail.com 
%              Institute of Applied Mechanics, 2012.08.20
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: Jiang & Shu; Efficient Implementation of Weighted ENO Schemes
% JCP. vol 126, 202-228 (1996)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Finite Difference Implementation of WENO5 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
L = 100*10^-9;              %device length in meters
num_cell = 10000;            % number of cells
p_initial =  10^27;        %initial hole density
p_mob = 2.0*10^-8;         %hole mobility

Va_min = 1000;             %volts
Va_max = 1000;    
increment = 1;       %for increasing V
num_V = floor((Va_max-Va_min)/increment)+1;

%Simulation parameters
w = .01;              %set up of weighting factor
tolerance = 10^-14;   %error tolerance       
constant_p_i = true;
fluxsplit = 3;        % {1} Godunov, {2} Global LF, {3} Local LF  Defines which flux splitting method to be used

%% Physical Constants
q =  1.60217646*10^-19;         %elementary charge, C
kb = 1.3806503D-23;              %Boltzmann const., J/k
T = 296.;                       %temperature
epsilon_0 =  8.85418782*10^-12; %F/m
epsilon = 3.8*epsilon_0;        %dielectric constant of P3HT:PCBM

Vt = (kb*T)/q;

%% Define our Flux function
fluxtype = 'linear';

%These @(w) define flux and dflux as function handles, allowing to call
%functions indirectly within arguments of calls to other functions as done
%in residual().
switch fluxtype
    case 'linear'
        c=1; flux = @(w) c*w;
        dflux = @(w) c*ones(size(w));
    case 'nonlinear' % Burgers'
        flux = @(w) w.^2/2;
        dflux = @(w) w;
end

%% Domain Discretization
a=0; b=L; x0=linspace(a,b,num_cell+1); dx=(b-a)/num_cell;   %x0 is positions array, +1 necessary b/c linspace needs (# of pts = # of cells +1)

% 2 more ghost pts on each side have to be added to the final domain.
x= [-2*dx,-dx,x0,b+dx,b+2*dx];   
nx = length(x);  

%% Initial Conditions
if(constant_p_i)
    for i = 3:nx-2
        p(i) = p_initial;
    end
else
    %linearly decreasing p: this doesn't work well
    p(3) = p_initial;
    for i = 3:nx-2      
        dp = p_initial/(num_cell+1);
        p(i+1) = p(i)-dp;
    end
end

%%Boundary conditions: define ghost points
p(1) = 0;
p(2) = 0;
p(nx-1) = 0;
p(nx) = 0;

    Va_cnt = 1;
for Va_cnt = 1:num_V
    
    %% Initial Conditions
    if(constant_p_i)
        for i = 3:nx-2
            p(i) = p_initial;
        end
    else
        %linearly decreasing p: this doesn't work well
        p(3) = p_initial;
        for i = 3:nx-2      
            dp = p_initial/(num_cell+1);
            p(i+1) = p(i)-dp;
        end
    end

    Va = Va_min+increment*(Va_cnt-1);    %increase Va
    %Va = Va_max-increment*(Va_cnt-1);    %decrease Va by increment in each iteration
    
    %for low Va_max (start point), use lower 1st w, for medium Va, use
    %lower w.
%     if(Va_max<200.)
%         if(Va_cnt==1)
%             w= 0.0001;
%         elseif(Va_max<30.)  %need to figureout what this value limit is
%             w = 0.0001;
%         else
%             w = 0.001;
%         end
%     elseif(Va<200.)
%         w = 0.001;
%     end
    
    %timing
    tic
    
    %% Setup the solver
    iter = 0;
    error_p =  1.0;
    
    %set up AV: these never change so do outside of loop
    %set up using sparse matrix (just store the lower diag, main diag,upper
    %diag in num_pts x 3 matrix AV_val.
    num_pts = nx-3-4+1;          %(=# number nodes inside device +1);  
    AV = zeros(num_pts,3);
    AV(:,1) = 1.;
    AV(:,3) = 1.;
    AV(:,2) = -2.;
    AV_val = spdiags(AV,-1:1,num_pts,num_pts);  %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
    
    %allocate matrices/arrays
    V = zeros(num_pts,1);
    Ap =  zeros(num_pts,3);
    bV = zeros(num_pts,1);
    bp = zeros(num_pts,1);
    
    %% Solver Loop                
    while error_p > tolerance
        
        %Poisson equation with tridiagonal solver
    
        % setup bV
        for i = 1: num_pts
            bV(i,1) = -dx^2*(q/epsilon)*p(i+3);    %bV(i) corresponds to 1st point inside device which is 4 so i+3
        end
        %for bndrys (AV and bV will be solved on range 4:nx-3)
        bV(1,1) = bV(1,1) - Va;
        bV(num_pts,1) = bV(num_pts,1) - 0;
        
        %call solver, solve for V
        V =  AV_val\bV;
        
        %make V with ghost points and bndry pts
        fullV = [0; 0; Va; V; 0; 0; 0];   %THIS SHOULD BE FORCED TO 0 AT RIGHT BOUNDARY! FIND THAT THIS GIVES SMOOTHER CURVE AND IS CORRECT
        %NOT SURE IF ALL LEFT POINTS SHOULD BE Va or should left 2 points
        %be 0: doesn't seem to make much  difference but probably b/c not

        fullV = fullV.';             %transpose to be able to put into residual
        
        %finite differences: for comparison
        for i = 3:nx-3
            dV(i) = (fullV(i+1)-fullV(i))/dx;
        end  
  
        %ALL OF THE BELOW ARE NECESSARY FOR THE RESULTS TO BE CORRECT
        %tested by removing some of them: then really has issues.
        dV(3) = (fullV(4)-fullV(3))/dx;   
        dV(4) = (fullV(5)-fullV(4))/dx; 
        dV(nx-2) = (fullV(nx-2)-fullV(nx-3))/dx;
        dV(nx-3) = (fullV(nx-3)-fullV(nx-4))/dx; 
      
        E = -dV;
        
        %BCs: set E's to equal the values that they are right inside the
        %device
        E(1) = E(3);
        E(2) = E(3);
        E(nx) = E(nx-2);
        E(nx-1) = E(nx-2);
        
        %% now solve eqn for p

        %finite differences for dp and dE
        for i = 3:nx-3
            dE(i) = (E(i+1)-E(i))/dx;
            dp(i) = (p(i+1)-p(i))/dx;
        end
        
        %Take care of boundaries:
        dE(3) = (E(4)-E(3))/dx;     %IF DON'T IMPOSE THIS BC, IT gives completely wrong results    
        dE(nx-2) = (E(nx-2)-E(nx-3))/dx;
        dE(nx-1) = 0;
        
        dp(3) = (p(4)-p(3))/dx;
        dp(nx-2) = (p(nx-2)-p(nx-3))/dx;
        dp(nx-1) = 0;
           
        old_p = p;
        
        %setup matrices
        Ap = zeros(num_pts,3);
        Ap(:,1) = 1.;
        Ap(:,3) = 1.;
    
        Cp = dx^2/Vt;
        for k = 1:num_pts
            Ap(k,2) = -2.-(dE(k+3))/Vt;    %dE(4) corresponds to 1st point within device which will correspond to 1st element in matrix/arrays
            bp(k) = Cp*E(k+3)*dp(k+3);
        end
        bp(1) = bp(1) - p(3);
        bp(num_pts) = bp(num_pts) - p(num_pts+1);
        Ap_val = spdiags(Ap,-1:1,num_pts,num_pts); %POTENTIAL ISSUE!! THE MAIN DIAGONAL ELEMENTS OF THIS MATRIX ARE 10^18 larger than upper and lower diag elements b/c of the dE/dx
        
        p_sol = Ap_val\bp;   
        fullp = [0;0;p_initial; p_sol;0; 0;0];
        newp = fullp.';  %transpose so matches with other matrices (horizontal array, 1 row).
        
        %weighting
        %p = newp;
        p = newp*w + old_p*(1-w);
        
        error_p = max(abs(p-old_p)/abs(old_p))
    
       iter =  iter+1    
       
       if(Va == Va_min) %only for last run
          E_solution(iter) = E(46);    % save E at random point for each iter for convergence analysis
       end
    end
    
    %Calculate drift diffusion Jp
    %update dp
    for i = 3:nx-3
         dp(i) = (p(i+1)-p(i))/dx;
    end
    dp(3) = (p(4)-p(3))/dx;
    dp(nx-2) = (p(nx-2)-p(nx-3))/dx;
    dp(nx-1) = 0;
    
    for i = 3:nx-2 
        Jp(i) = q*p_mob*p(i)*E(i)- q*p_mob*Vt*dp(i);
    end
    
    %Setup for JV curve
    V_values(Va_cnt) = Va;
    Jp_final(Va_cnt) = Jp(nx-3);  %just pick Jp at the right side
       
    %Save data
    str = sprintf('%.2f',Va);
    filename = [str 'V.txt'] 
    fid = fopen(fullfile('C:\Users\Tim\Documents\Duxbury group research\WENO\results_output\',filename),'w');   %w: Open or create new file for writing
    %fullfile allows to make filename from parts
    for i = 3:nx-2
        fprintf(fid,'%.8e %.8e %.8e %.8e\r\n ',x(i), p(i), E(i), Jp(i)); 
        %f means scientific notation, \r\n: start new line for each value
        %for each new desired column, put its own % sign    
    end
    fclose(fid);
    
    toc
    
end

%sanity check: calculate V by integrating E
% V_final(Va_cnt) = 0;
%     for i = 3:nx-2
%         V_final(Va_cnt) = V_final(Va_cnt) + E(i)*dx;
%     end  
%     V_final

%% Final Plots

str = sprintf('%.2g', Va);

 h1 = plot(x(3:num_cell+3),p(3:num_cell+3));
 hold on
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);
 
 figure;
 h2 = plot(x(3:num_cell+3),E(3:num_cell+3));
 hold on
 %plot(x(3:num_cell),E_theory1(3:num_cell));
 %plot(x(3:num_cell),E_theory2(3:num_cell));
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);
 
 
 figure;
 h3 = plot(x(3:num_cell+3),Jp(3:num_cell+3));
 hold on
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
 
 %JV curve
 figure
 h4 = plot(V_values,Jp_final);
 hold on
 %plot(V_values, Jp_theory);
 xlabel('Voltage (V)','interpreter','latex','FontSize',14);
 ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
 
 %convergence analysis
 iterations = 1:iter;
 figure
 h5 = plot(iterations, E_solution);
 title(['E convergence', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Iterations','interpreter','latex','FontSize',14);
 ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);
 
 hold off
