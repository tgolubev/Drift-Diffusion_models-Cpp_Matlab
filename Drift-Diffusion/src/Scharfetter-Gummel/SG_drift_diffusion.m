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
num_cell = 100;            % number of cells
p_initial =  10^27;        %initial hole density
p_mob = 2.0*10^-8;         %hole mobility

Va_min = 1;             %volts
Va_max = 1;    
increment = 1;       %for increasing V
num_V = floor((Va_max-Va_min)/increment)+1;

%Simulation parameters
w = 0.5;              %set up of weighting factor
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
a=0; b=L; x=linspace(a,b,num_cell+1); dx=(b-a)/num_cell;   %x is positions array, +1 necessary b/c linspace needs (# of pts = # of cells +1)  
nx = length(x);  

%% Initial Conditions
fullp = zeros(nx,1);     %make it vertical column vector for convienicne with matrix operations
if(constant_p_i)
    for i = 1:nx
        fullp(i,1) = p_initial;
    end
else

    %linearly decreasing p: this doesn't work well
%     fullp(3) = p_initial;
%     for i = 1:nx      
%         dp = p_initial/(num_cell+1);
%         p(i+1) = p(i)-dp;
%     end
end

    Va_cnt = 1;
for Va_cnt = 1:num_V
    
%     %% Initial Conditions
%     if(constant_p_i)
%         for i = 1:nx
%             p(i) = p_initial;
%         end
%     else
%         %linearly decreasing p: this doesn't work well
%         p(3) = p_initial;
%         for i = 1:nx      
%             dp = p_initial/(num_cell+1);
%             p(i+1) = p(i)-dp;
%         end
%     end

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
        
        %Poisson equation with tridiagonal solver
    
        % setup bV
        for i = 1: num_elements
            bV(i,1) = -dx^2*(q/epsilon)*fullp(i+1);   
        end
        %for bndrys 
        bV(1,1) = bV(1,1) - Va;
        bV(num_elements,1) = bV(num_elements,1) - 0;
        
        %call solver, solve for V
        V =  AV_val\bV;
     
        %make V with bndry pts
        fullV = [Va; V; 0];   %THIS SHOULD BE FORCED TO 0 AT RIGHT BOUNDARY! 

        fullV = fullV.';             %transpose to be able to put into residual
        
%------------------------------------------------------------------------------------------------        
       %% now solve eqn for p           
        %B = BernoulliFnc(nx, fullV, Vt);
        
        fullp(nx) = fullp(nx-1);
       
        %Ap_val = SetAp_val(num_elements, B, fullV,p,Vt);      
 for i = 1:nx-1               %so dV(100) = dV(101: is at x=L) - dV(100, at x= l-dx).
    dV(i) = fullV(i+1)-fullV(i);
    B(1,i) = dV(i)/(exp(dV(i))-1.0);    %B(+dV)
    B(2,i) = B(1)*exp(dV(i));          %B(-dV)
 end
        

 for i=1:num_elements
    Ap(i,1) = B(1,i);     %lower diagonal
    Ap(i,2) = -(B(2,i) + B(1,i+1));  %main diagonal
    Ap(i,3) = B(2,i+1);     %upper diagonal   CHECK MAKE SURE THIS IS RIGHT!
 end

 
Ap_val = spdiags(Ap,-1:1,num_elements,num_elements); %A = spdiags(B,d,m,n) creates an m-by-


        %enforce boundary conditions through bp
        bp(1) = bp(1) - B(1,1)*fullp(1);
        bp(num_elements) = bp(num_elements) - B(2,nx-1)*fullp(nx);           %NOTE I'M USING DIFFERENT NOTATION THAN DD-BI CODE!!: B's are defined from the left side!

        p_sol = Ap_val\bp;   

       error_p = max(abs(p_sol-fullp(2:nx-1))/abs(fullp(2:nx-1)));
        
        %weighting
        oldp = fullp(2:nx-1);
        p = p_sol*w + oldp*(1-w);      %fullp is soln. from  previous iter
  
        
        fullp = [p_initial; p];  %leave right side unforced
        %fullp = fullp.';  %transpose so matches with other matrices (horizontal array, 1 row).

        iter =  iter+1    
       
       if(Va == Va_min) %only for last run
          p_solution(iter) = fullp(46);    % save E at random point for each iter for convergence analysis
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
    for i = 1:nx-1
        fprintf(fid,'%.8e %.8e\r\n ',x(i), fullp(i)); 
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

plotp = fullp.';   %transpose to be able to plot
 h1 = plot(x(1:nx-1),plotp);
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