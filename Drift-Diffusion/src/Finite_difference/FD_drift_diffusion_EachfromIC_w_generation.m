%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving Poisson + Drift Diffusion eqns (for holes) with finite 
%           differences. 
%
%NOTE: RIGHT NOW THIS CAN'T CONVERGE!
%
%                  Coded by Timofey Golubev (2017.08.11)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
L = 100*10^-9;              %device length in meters
num_cell = 100;            % number of cells
p_initial =  10^23;        %initial hole density
p_mob = 2.0*10^-8;         %hole mobility

Va_min = 1;             %volts
Va_max = 1;    
increment = 1;       %for increasing V
num_V = floor((Va_max-Va_min)/increment)+1;

%Simulation parameters
U = 0;%10^25;           %net hole generation rate
w = 0.01;              %set up of weighting factor
tolerance = 10^-14;   %error tolerance       
constant_p_i = true;

%% Physical Constants
q =  1.60217646*10^-19;         %elementary charge, C
kb = 1.3806503D-23;              %Boltzmann const., J/k
T = 296.;                       %temperature
epsilon_0 =  8.85418782*10^-12; %F/m
epsilon = 3.8*epsilon_0;        %dielectric constant of P3HT:PCBM

Vt = (kb*T)/q;

%% Domain Discretization
a=0; b=L; x=linspace(a,b,num_cell+1); dx=(b-a)/num_cell;   %x0 is positions array, +1 necessary b/c linspace needs (# of pts = # of cells +1) 
nx = length(x);  

%% Initial Conditions
if(constant_p_i)
    for i = 1:nx
        p(i) = p_initial;
    end
else
    %linearly decreasing p: this doesn't work well
    p(1) = p_initial;
    for i = 1:nx      
        dp = p_initial/(num_cell+1);
        p(i+1) = p(i)-dp;
    end
end

  %redefine p's to be only those inside device
    p = p(2:num_cell);
    Va_cnt = 1;
for Va_cnt = 1:num_V
    
    Va = Va_min+increment*(Va_cnt-1);    %increase Va
    %Va = Va_max-increment*(Va_cnt-1);    %decrease Va by increment in each iteration
    
    %timing
    tic
    
    %% Setup the solver
    iter = 0;
    error_p =  1.0;
    
    %set up AV: these never change so do outside of loop
    %set up using sparse matrix (just store the lower diag, main diag,upper
    %diag in num_pts x 3 matrix AV_val.
    
    %Here don't have to worry about the exact setup of colunms for sparse
    %matrix b/c all values in each diagonal are same anyway.
    AV = zeros(nx-2,3);
    AV(:,1) = 1.;
    AV(:,3) = 1.;
    AV(:,2) = -2.;
    AV_val = spdiags(AV,-1:1,nx-2,nx-2);  %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
    
    %allocate matrices/arrays
    V = zeros(nx-2,1);
    Ap =  zeros(nx-2,3);
    bV = zeros(nx-2,1);
    bp = zeros(nx-2,1);
    
    %% Solver Loop 
    num_elements = nx-2;
    CV = dx^2*q/epsilon;
    while error_p > tolerance
        
        %Poisson equation with tridiagonal solver
        % setup bV
        for i = 1:num_elements
            bV(i,1) = -CV*p(i);    
        end
        %for bndrys (AV and bV will be solved on range 4:nx-3)
        bV(1,1) = bV(1,1) - Va;
        bV(num_elements,1) = bV(num_elements,1) - 0;
        
        %call solver, solve for V
        V =  AV_val\bV;
     
        %make V with bndry pts
        fullV = [Va; V; 0];   %THIS SHOULD BE FORCED TO 0 AT RIGHT BOUNDARY! FIND THAT THIS GIVES SMOOTHER CURVE AND IS CORRECT
        fullV = fullV.';             %transpose to be able to put into residual
         
        %finite differences: for comparison
        for i = 1:nx-1
            dV(i) = (fullV(i+1)-fullV(i))/dx;
        end  
 
        E = -dV;
        
        %BCs: set right side E to equal the values right inside the bndry
        E(nx) =  E(nx-1);
        
        %% now solve eqn for p
        
        fullp = [p_initial, p, 0];  %add bndry values (right side bndry p =0)

        %finite differences for dp and dE
        for i = 1:num_elements
            dE(i) = (E(i+1)-E(i))/dx;
            dp(i) = (fullp(i+1)-fullp(i))/dx;
        end
           
        old_p = p;
        
        %setup matrices
        %Don't have to worry about right lengths of lower and upper diags
        %here, since all elements are 1.
        Ap = zeros(num_elements,3);
        Ap(:,1) = 1.;
        Ap(:,3) = 1.;
    
        Cp = dx^2/Vt;
        for k = 1:num_elements
            Ap(k,2) = -2.-(dE(k))/Vt;    %dE(4) corresponds to 1st point within device which will correspond to 1st element in matrix/arrays
            bp(k) = Cp*E(k+1)*dp(k);
        end
        bp(1) = bp(1) - fullp(1);
        
        %add generation in approx middle
        bp(floor(num_cell/2.)) = bp(floor(num_cell/2.)) + U/(Vt*p_mob);
        
        bp(num_elements) = bp(num_elements) - fullp(num_elements+1);   %this fullp value  at right  side = 0
        Ap_val = spdiags(Ap,-1:1,num_elements,num_elements); %POTENTIAL ISSUE!! THE MAIN DIAGONAL ELEMENTS OF THIS MATRIX ARE 10^18 larger than upper and lower diag elements b/c of the dE/dx
        
        p_sol = Ap_val\bp;   

        newp = p_sol.';  %transpose so matches with other matrices (horizontal array, 1 row).
        
        error_p = max(abs(newp-old_p)/abs(old_p))   %should calculate error before weighting
        
        %weighting
        %p = newp;
        p = newp*w + old_p*(1-w);
 
       iter =  iter+1    
       
       if(Va == Va_min) %only for last run
          E_solution(iter) = E(46);    % save E at random point for each iter for convergence analysis
       end
    end
    
    %Calculate drift diffusion Jp
    %update dp
    for i = 1:num_elements
         dp(i) = (p(i+1)-p(i))/dx;
    end
    
    for i = 1:num_elements
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
    fullp = [p_initial, p];  %add left side value
    for i = 1:nx-2
        fprintf(fid,'%.8e %.8e\r\n ',x(i+1), fullp(i));   %x(i+1) b/c not printing bndry p's
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

 h1 = plot(x(1:nx-1),fullp);
 hold on
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);
 
 %Plot potential (fullV)
 figure
 plot(x, fullV)
 hold on
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Electric Potential (V)'},'interpreter','latex','FontSize',14);
 
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
 h5 = plot(iterations, E_solution);
 title(['E convergence', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Iterations','interpreter','latex','FontSize',14);
 ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);
 
 hold off
