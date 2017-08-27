%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving Poisson + Drift Diffusion eqns (for holes) using
%                   Scharfetter-Gummel discretization
%
%                  Coded by Timofey Golubev (2017.08.16)
%             NOTE: i=1 corresponds to x=0, i=nx to x=L

%RECALL THAT fullV = V/Vt: is scaled!. Same with p = p/N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;   %NOTE: clear improves performance over clear all, and still clears all variables.

%% Parameters
L = 100*10^-8;              %device length in meters   %NOTE IIF USE 100*10^-8 GET CORRECT MOTT -GURNEY LOOKING BEHAVIOR. IF 100*10^-9: get weird almost flat p.
                            % b/c of V falling off too fast.
num_cell = 100;            % number of cells     %IF USE 100, get kink in left side of E_theory. At 1000 cells, kink is gone.
p_initial =  10^23;        %initial hole density   %NOTE: WORKS FOR UP TO 10^23, BEYOND THAT, HAVE ISSUES
p_mob = 2.0*10^-8;         %hole mobility

U = 0;%10^28;                       %net carrier generation rate at interface (in middle): NOTE: BOTH MINUS AND PLUS WORK! (minus up to  -10^30).

N = 10^23.;   %scaling factor helps CV be on order of 1 but makes Cp be very small --> not sure how useful it is, unless figure out how to better scale Cp also..

p_initial = p_initial/N;

Va_min = 35;             %volts       NOTE: I;M NOW RELAXING TOLERANCE FOR 1ST ITER, SO CAN'T TRUST VA_MIN RESULT: RUN FOR AT LEAST 2 Va values!
Va_max = 35.1;    
increment = 0.1;       %for increasing V
num_V = floor((Va_max-Va_min)/increment)+1;

%Simulation parameters
w = 0.05;              %set up of weighting factor     %IT seems to WORK WITH w = 0.5 without issues for 100nodes:  WHEN USING MORE NODES, NEED LOWER w: like 0.01
tolerance = 10^-13;   %error tolerance       
constant_p_i = true;   


%% Physical Constants
q =  1.60217646*10^-19;         %elementary charge, C
kb = 1.3806503*10^-23;              %Boltzmann const., J/k
T = 296.;                       %temperature
epsilon_0 =  8.85418782*10^-12; %F/m
epsilon = 3.8*epsilon_0;        %dielectric constant of P3HT:PCBM

Vt = (kb*T)/q;

%% Domain Discretization
a=0; b=L; x=linspace(a,b,num_cell+1); dx=(b-a)/num_cell;   %x is positions array, +1 necessary b/c linspace needs (# of pts = # of cells +1)  
nx = length(x);  

%% Initial Conditions
if(constant_p_i)
    for i = 1:nx
        p(i) = 10^-20; % make p initial very small--> this corresponds to how it's done in ddbi code    ACUTALLY FIND THAT IT DOESN'T MATTER WHAT THIS SET TO--> CAN SET TO 1 AND STILL WORKS.
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
    
    if(Va_cnt ==1) 
        tolerance = tolerance*1000;         %relax tolerance for 1st convergence
    end
    if(Va_cnt ==2)
        tolerance = tolerance/1000;       %reset tolerance back
    end

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
    
   %Here don't have to worry about the exactl setup of colunms for sparse
   %matrix b/c all values in each diagonal are same anyway.
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
    Cp = dx^2/(Vt*N*p_mob);   %note: I divided the p_mob out of the matrix
    CV = N*dx^2*q/(epsilon*Vt);
    while error_p > tolerance
         
        %Poisson equation with tridiagonal solver
    
        % setup bV
        for i = 1: num_elements
            bV(i,1) = -CV*p(i);   
        end
        %for bndrys 
        bV(1,1) = bV(1,1) - Va/Vt;         %must scale this too, since all others are scaled
        bV(num_elements,1) = bV(num_elements,1) - 0;
        
        %call solver, solve for V
        %spparms('spumoni',2)     %THIS ALLOWS TO SEE THE DETAILS OF THE SOLVER CONDITIONS!
        V =  AV_val\bV;
        %SO THIS ACTUALLY IS SOLVING FOR PSI PRIME: A SCALED PSI! --> so
        %later in Bernoulli: don't need to devide by VT again: b/c this
        %here is already scaled!
        
        %make V with bndry pts
        fullV = [Va/Vt; V; 0];   %THIS SHOULD BE FORCED TO 0 AT RIGHT BOUNDARY! 

        fullV = fullV.';             %transpose to be able to put into residual
  
     %NOTE: IT IS WRONG TO SCALE fullV so don't do that!
     
%------------------------------------------------------------------------------------------------        
       %% now solve eqn for p  
        %B = BernoulliFnc(nx, fullV, Vt);
 
        %Ap_val = SetAp_val(num_cell, B, fullV,p,Vt);      
        for i = 2:num_cell+1               %so dV(100) = dV(101: is at x=L) - dV(100, at x= l-dx).
            dV(i) = fullV(i)-fullV(i-1);     %these fullV's ARE ACTUALLY PSI PRIME; ALREADY = V/Vt, when solved poisson.
        end
        
        for i = 2:num_cell+1
            B(1,i) = dV(i)/(exp(dV(i))-1.0);    %B(+dV)
            %B(2,i) = -dV(i)/(exp(-dV(i))-1.0);   %THIS IS EQUIVALENT  TO OTHE
            %BELOW EXPERSSION
            B(2,i) = B(1,i)*exp(dV(i));          %B(-dV)
        end
        
        %BE CAREFUL!: the setup of the  sparse matrix elements is
        %non-trivial: last entry  of Ap(:,1) 1st entry of
        %Ap(:,3) are unused b/c off diagonals have 1 less element than main
        %diagonal!
        for i=1:num_elements     %I verified that ordering of columns is correct!
            Ap(i,2) = -(B(2,i+1) + B(1,i+1+1));     %I HAVE VERIFIED THAT ALL THE dV's properly  match up with  indices! All the extra +1's are b/c B starts at i=2 (1st pt inside device).
        end
        for i = 1:num_elements-1
            Ap(i,1) = B(1,i+1+1);        %THERE WAS INDEX MISTAKE HERE!  HAD B(1,i+1) but should be i+1+1 b/c have B(dVi)*p_i-1  (see fortran reindexed version)
        end
        Ap(num_elements, 1) = 0;   %last element is unused
        for i = 2:num_elements
            Ap(i,3) = B(2,i+1+1);   %WAS INDEX MISTAKE HERE ALSO
        end
        Ap(1,3) = 0;        %1st element of Ap(:,3) is unused
        
        
        Ap_val = spdiags(Ap,-1:1,num_elements,num_elements); %A = spdiags(B,d,m,n) creates an m-by-
        old_p = p;
        
        %enforce boundary conditions through bp
        bp(1) = -B(1,2)*p_initial;       %NOTE: scaled p_initial = 1 --> this is just like in fortran version
        bp(num_elements) = 0;%-B(1,nx-1)*p(nx);       %ENFORCE RIGHT SIDE P IS = 0 that set  %NOTE I'M USING DIFFERENT NOTATION THAN DD-BI CODE!!: B's are defined from the left side!
        %dmaybe here need to include THE -B for the other boundary c
        %ondition: subtract the right side p value.
        
        %introduce a net generation rate somewhere in the middle
        bp(ceil(num_cell/2.)) = -Cp*U;
            
        p_sol = Ap_val\bp;
   
        newp = p_sol.';    %tranpsose
        error_p = max(abs(newp-old_p)/abs(old_p))  %ERROR SHOULD BE CALCULATED BEFORE WEIGHTING
        
        %weighting
        p = newp*w + old_p*(1.-w);
   
       iter =  iter+1    
  
       if(Va == Va_max) %only for last run
          p_solution(iter) = p(20);    % save E at random point for each iter for convergence analysis
       end
    end
    %Define fullp
    fullp = [p_initial, p, 10^-20];  %add bndry values: right side set to very small value --> same as did in Fortran version
    
    %     %Calculate drift diffusion Jp
    % Use the SG definition
    for i = 1:num_cell-1        %verified that this is right 
        Jp_temp(i) = -(q*p_mob*Vt*N/dx)*(fullp(i+1)*B(2,i+1)-fullp(i)*B(1,i+1));         %need an N b/c my p's are scaled by N. Extra +1's are b/c B starts at i=2
    end
    
    for i =  2:num_cell
        Jp(i) = Jp_temp(i-1);     %define Jp as on the right side (i.e. i+1/2 goes to i+1).
    end
    
    %Setup for JV curve
    V_values(Va_cnt) = Va;
    Jp_final(Va_cnt) = Jp(nx-3);  %just pick Jp at the right side
    
    %numerical E result calculation
    for i = 2:num_cell+1
        E(i) = -Vt*(fullV(i) - fullV(i-1))/dx;  %*Vt to rescale  back to normal units: RECALL THAT I DEFINED dV as just V(i+1) - V(i) and here we need dV/dx!!
    end
    
    %Analytic Result Calculation
    for i=2:num_cell
        %E_theory1(i) = sqrt(2*x(i)*abs(Jp(nx-3))/(epsilon*p_mob));
        E_theory2(i)= sqrt(2*x(i)*abs(Jp(i))/(epsilon*p_mob));        %THIS IS MORE CORRECT WAY, SINCE USE THE Jp at each point
    end
    
 
    %Save data
    str = sprintf('%.2f',Va);
    filename = [str 'V.txt'] 
    fid = fopen(fullfile('C:\Users\Tim\Documents\Duxbury group research\WENO\results_output\',filename),'w');   %w: Open or create new file for writing
    %fullfile allows to make filename from parts
  
    for i = 2:num_cell
        fprintf(fid,'%.8e %.8e %.8e %.8e %.8e %.8e\r\n ',x(i), Vt*fullV(i), fullp(i), E(i), E_theory2(i), Jp(i) );   %x(i+1) b/c not printing bndry p's
        %f means scientific notation, \r\n: start new line for each value
        %for each new desired column, put its own % sign    
    end
    fclose(fid);
    
    toc
  
end

%output to screen for testing
p(10);
p(11);

%sanity check: calculate V by integrating E
% V_final(Va_cnt) = 0;
%     for i = 3:nx-2
%         V_final(Va_cnt) = V_final(Va_cnt) + E(i)*dx;
%     end  
%     V_final

%% Final Plots: done for the last Va

%Analytic Result Calculation
for i=2:num_cell
   %E_theory1(i) = sqrt(2*x(i)*abs(Jp(nx-3))/(epsilon*p_mob));
   E_theory2(i)= sqrt(2*x(i)*abs(Jp(i))/(epsilon*p_mob));        %THIS IS MORE CORRECT WAY, SINCE USE THE Jp at each point

end

str = sprintf('%.2g', Va);

 h1 = plot(x,fullp);
 hold on
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);
 axis([-inf inf 0 inf]);
 
 %Plot E
 figure
 E = -(dV/dx)*Vt;  %*Vt to rescale  back to normal units: RECALL THAT I DEFINED dV as just V(i+1) - V(i) and here we need dV/dx!!
 plot(x(2:num_cell), E(2:num_cell))      %We don't plot left bndry pt., b/c E there is not calculated. dV starts at i=2.
 hold on
 plot(x(2:num_cell),E_theory2(2:num_cell));
 %plot(x(2:nx-1),E_theory1);
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);
 if -dV(1) < 0.0
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

 
 
 figure;
 h3 = plot(x(2:num_cell),Jp(2:num_cell));
 hold on
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
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
