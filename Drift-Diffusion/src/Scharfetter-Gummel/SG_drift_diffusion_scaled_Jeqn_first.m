%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving Poisson + Drift Diffusion eqns (for holes) using
%                   Scharfetter-Gummel discretization
%
%                  Coded by Timofey Golubev (2017.08.11)
%             NOTE: i=1 corresponds to x=0, i=nx to x=L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
L = 100*10^-9;              %device length in meters
num_cell = 20000;            % number of cells
p_initial =  10^27;        %initial hole density
p_mob = 2.0*10^-8;         %hole mobility

N = 10^27;    %scaling factor for p

p_initial = p_initial/N;

Va_min = 1;             %volts
Va_max = 1;    
increment = 0.01;       %for increasing V
num_V = floor((Va_max-Va_min)/increment)+1;

%Simulation parameters
w = 0.001;              %set up of weighting factor
tolerance = 10^-14;   %error tolerance       
constant_p_i = true;   % NEED LINEARLY DECREASING P INITIALLY TO HAVE BETTER INITIAL V RESULT--> better stability


%% Physical Constants
q =  1.60217646*10^-19;         %elementary charge, C
kb = 1.3806503D-23;              %Boltzmann const., J/k
T = 296.;                       %temperature
epsilon_0 =  8.85418782*10^-12; %F/m
epsilon = 3.8*epsilon_0;        %dielectric constant of P3HT:PCBM

Vt = (kb*T)/q;

%% Domain Discretization
a=0; b=L; x=linspace(a,b,num_cell+1); dx=(b-a)/num_cell;   %x is positions array, +1 necessary b/c linspace needs (# of pts = # of cells +1)  
nx = length(x);  

%% Initial Conditions


dv = Va_min/(Vt*(num_cell+1));
fullV(1) = Va_min/Vt;
for i = 1:nx-1
    fullV(i+1) = fullV(i) - dv;
    p(i) = p_initial;                    %make p_initial just constantt throughout
end


    Va_cnt = 1;
for Va_cnt = 1:num_V
    

    
  %redefine p's to be only those inside device
    p = p(2:num_cell);

    Va = Va_min+increment*(Va_cnt-1);    %increase Va
    %Va = Va_max-increment*(Va_cnt-1);    %decrease Va by increment in each iteration
    
    %timing
    tic
    
    %% Setup the solver
    iter = 1;
    error_p =  1.0;
    
    %set up AV: these never change so do outside of loop
    %set up using sparse matrix (just store the lower diag, main diag,upper
    %diag in num_pts x 3 matrix AV_val.
    
    AV = zeros(nx-2,3);
    AV(:,1) = 1.;
    AV(:,3) = 1.;
    AV(:,2) = -2.;
    AV_val = spdiags(AV,-1:1,nx-2,nx-2);  %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
    
    %allocate matrices/arrays
    V = zeros(nx-2,1);
    bV = zeros(nx-2,1);
    
    B = zeros(2, nx-1);
    bp = zeros(nx-2,1);
    
    %% Solver Loop        
    num_elements = nx-2;
    while error_p > tolerance
        
        
        %% Solve eqn for p first        
        %B = BernoulliFnc(nx, fullV, Vt);
       
        %Ap_val = SetAp_val(num_cell, B, fullV,p,Vt);      
 for i = 1:nx-1               %so dV(100) = dV(101: is at x=L) - dV(100, at x= l-dx).
    dV(i) = fullV(i+1)-fullV(i);     %these fullV's ARE ACTUALLY PSI PRIME; ALREADY = V/Vt, when solved poisson.
 end
 %add injection step
 dV(1) = dV(1) + 1./Vt;    %the # = phi_a in ddbi code
    
 for i = 1:nx-1
    B(1,i) = dV(i)/(exp(dV(i))-1.0);    %B(+dV)
    B(2,i) = B(1)*exp(dV(i));          %B(-dV)
 end


 for i=1:num_elements     %I verified that ordering of columns is correct!
    Ap(i,1) = B(1,i);    
    Ap(i,2) = (B(2,i) + B(1,i+1));  
    Ap(i,3) = B(2,i+1);    
 end


Ap_val = spdiags(Ap,-1:1,num_elements,num_elements); %A = spdiags(B,d,m,n) creates an m-by-
         old_p = p;
        
        %enforce boundary conditions through bp
        bp(1) = B(1,1)*p_initial;
        bp(num_elements) = 0;       %ENFORCE RIGHT SIDE P IS 0    %NOTE I'M USING DIFFERENT NOTATION THAN DD-BI CODE!!: B's are defined from the left side!
   
        p_sol = Ap_val\bp;   
        
        %fullp = [p_initial; p_sol;0];
        %newp = fullp.';  %transpose so matches with other matrices (horizontal array, 1 row).
        
        %weighting
        newp = p_sol.';    %tranpsose
        p = newp*w + old_p*(1.-w);
      
        error_p = max(abs(p-old_p)/abs(old_p))
        
      
        
        %Poisson equation with tridiagonal solver
    
        % setup bV
        for i = 1: num_elements
            bV(i,1) = -N*dx^2*(q/(epsilon*Vt))*p(i);   
        end
        %for bndrys 
        bV(1,1) = bV(1,1) - Va/Vt;         %must scale this too, since all others are scaled
        bV(num_elements,1) = bV(num_elements,1) - 0;
        
        %call solver, solve for V
        V =  AV_val\bV;
        %SO THIS ACTUALLY IS SOLVING FOR PSI PRIME: A SCALED PSI! --> so
        %later in Bernoulli: don't need to devide by VT again: b/c this
        %here is already scaled!
        
        %make V with bndry pts
        fullV = [Va/Vt; V; 0];   %THIS SHOULD BE FORCED TO 0 AT RIGHT BOUNDARY! 

        fullV = fullV.';             %transpose to be able to put into residual 

       iter =  iter+1    
  
       if(Va == Va_max) %only for last run
          p_solution(iter) = p(20);    % save E at random point for each iter for convergence analysis
       end
    end
    
%     %Calculate drift diffusion Jp
%     %update dp
%     for i = 3:nx-3
%          dp(i) = (p(i+1)-p(i))/dx;
%     end
%     dp(3) = (p(4)-p(3))/dx;
%     dp(nx-2) = (p(nx-2)-p(nx-3))/dx;
%     dp(nx-1) = 0;
    
%     for i = 3:nx-2 
%         Jp(i) = q*p_mob*p(i)*dV(i)- q*p_mob*Vt*dp(i);
%     end
%     
%     %Setup for JV curve
%     V_values(Va_cnt) = Va;
%     Jp_final(Va_cnt) = Jp(nx-3);  %just pick Jp at the right side
       
    %Save data
    str = sprintf('%.2f',Va);
    filename = [str 'V.txt'] 
    fid = fopen(fullfile('C:\Users\Tim\Documents\Duxbury group research\WENO\results_output\',filename),'w');   %w: Open or create new file for writing
    %fullfile allows to make filename from parts
    for i = 3:nx-2
        fprintf(fid,'%.8e %.8e\r\n ',x(i), p(i)); 
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

 h1 = plot(x(2:nx-1),p);
 hold on
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);
 
%  figure;
%  h2 = plot(x(3:num_cell+3),E(3:num_cell+3));
%  hold on
%  %plot(x(3:num_cell),E_theory1(3:num_cell));
%  %plot(x(3:num_cell),E_theory2(3:num_cell));
%  title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
%  xlabel('Position ($m$)','interpreter','latex','FontSize',14);
%  ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);
%  
%  
%  figure;
%  h3 = plot(x(3:num_cell+3),Jp(3:num_cell+3));
%  hold on
%  title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
%  xlabel('Position ($m$)','interpreter','latex','FontSize',14);
%  ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
%  
%  %JV curve
%  figure
%  h4 = plot(V_values,Jp_final);
%  hold on
%  %plot(V_values, Jp_theory);
%  xlabel('Voltage (V)','interpreter','latex','FontSize',14);
%  ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
 
 %convergence analysis
 iterations = 1:iter;
 figure
 h5 = plot(iterations, p_solution);
 title(['p convergence', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Iterations','interpreter','latex','FontSize',14);
 ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);
 
 hold off