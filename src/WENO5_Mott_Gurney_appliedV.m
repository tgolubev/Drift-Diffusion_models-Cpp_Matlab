%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving Poisson + Mott-Gurney limit system (for holes) with 5th order
%             Weighted Essentially Non-Oscilatory (WENO5)
%
%                E*dp/dx + p*dE/dx = 0  for x in [0,L]
%                  
%            
% Modified by Timofey Golubev (2017.08.06) based on original 1D wave eqn
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
num_cell = 100;            % number of cells
p_initial =  10^27;        %initial hole density
p_mob = 2.0*10^-8;         %hole mobility


Va_min = 200;             %volts
Va_max = 200;    
increment = 1.;        %for increasing V
%Ea_min = Va_min/L;         %V/m
%Ea_max = Va_max/L;         %maximum applied E
%increment = V_increment/L; %for increasing E

%Simulation parameters
w = .001;              %weighting factor
tolerance = 10^-14;   %error tolerance       
constant_p_i = true;
fluxsplit = 3;        % {1} Godunov, {2} Global LF, {3} Local LF  Defines which flux splitting method to be used

%% Physical Constants
q =  1.60217646*10^-19;         %elementary charge, C
epsilon_0 =  8.85418782*10^-12; %F/m
epsilon = 3.8*epsilon_0;        %dielectric constant of P3HT:PCBM

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
%matlab indices can't start with 0 so p(x=0) corresponds to i=1
% for i = 1:3
%     p(i) = p_initial;
% end


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

    Va_cnt = 0;
for Va = Va_min:increment:Va_max
    Va_cnt = Va_cnt +1;
  
%timing
tic
    %% Solver Loop
    iter = 0;
    error_p =  1.0;
    
    %set up AV_val: these never change so do outside of loop
    %set up using sparse matrix (just store the lower diag, main diag,upper
    %diag in num_pts x 3 matrix AV_val.
    num_pts = nx-3-4+1; %(=# number nodes inside device +1);  %setup the matrix
    A = zeros(num_pts,3);
    A(:,1) = 1.;
    A(:,3) = 1.;
    A(:,2) = -2.;
    
    %setup  sparse matrix
    AV = spdiags(A,-1:1,num_pts,num_pts);  %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
    
    %allocate matrices/arrays
    V = zeros(num_pts,1);
    bV = zeros(num_pts,1);
    
    %This sets up an entire matrix: is inefficient
%     num_pts = nx-3-4+1; %(=# nodes inside bndrys +1) ;
%     %define indices as integers so i==j works
%     %i = int8;
%     %j = int8;
%     for i = 1:num_pts
%         for j  = 1:num_pts
%             if i==j
%                 AV_val(i,j) = -2.;
%             else if abs(i-j) == 1
%                 AV_val(i,j) = 1;
%             else
%                 AV_val(i,j) = 0;
%             end
%             end
%         end
%     end        
    %AV_val
                    
    while error_p > tolerance
        
        %Poisson equation with tridiagonal solver
    
        % setup bV
        for i = 1: num_pts
            bV(i,1) = -dx^2*(q/epsilon)*p(i+3);    %bV(i) corresponds to 1st point inside device which is 4 so i+3
        end
        %for bndrys (AV and bV will be solved on range 4:nx-3)
        bV(1,1) = bV(1,1) - Va;%V(3);
        bV(num_pts,1) = bV(num_pts,1) - 0;%V(nx-2);
        
        %call solver, solve for V
        V =  AV\bV;
        
        %make V with ghost points and bndry pts
        fullV = [0; 0; Va; V; 0; 0; 0];   %THIS SHOULD BE FORCED TO 0 AT RIGHT BOUNDARY! FIND THAT THIS GIVES SMOOTHER CURVE AND IS CORRECT
        %IF FORCE LEFT ghost POINTS TO Va, really  messes it up.
        
        %USING 0 FOR LEFT GHOST POINTS SEEMS TO WORK BETTER THAN USING Va
        
        
        %dV = weno approx for dV/dx
        %NEED TO TRANSPOSE FULLV BEFORE PUTTING INTO WENO ALGO!!: A.' is
        %synthax for transpose of A
        fullV = fullV.';
        dV = residual(fullV,flux,dflux,dx,nx,fluxsplit);     %this calculates the entire dE array
        %ALL OF THE BELOW ARE NECESSARY TO ENSURE GOOD CONVERGENCE: I
        %tested by removing some of them: then really has issues.
        dV(3) = (fullV(4)-fullV(3))/dx;    %taking care of boundary issue
        dV(4) = (fullV(5)-fullV(4))/dx; 
        dV(nx-2) = (fullV(nx-2)-fullV(nx-3))/dx;
        dV(nx-3) = (fullV(nx-3)-fullV(nx-4))/dx; %this is necessary to ensure smooth monotonous dV
       
        %below condition doensn't seem to do anything
        %dV(5) = (fullV(6)-fullV(5))/dx;  
   
        E = -dV;
        
        %BCs: set E's to equal the values that they are right  inside the
        %device
        E(1) = E(3);
        E(2) = E(3);
        E(nx) = E(nx-2);
        E(nx-1) = E(nx-2);
        
        %now solve eqn for p
        dE = residual(E,flux,dflux,dx,nx,fluxsplit);
        %RIGHT NOW HAVING ISSUE THAT THE dE(3) is too large: having
        %issue doing the boundary derivative estimate properly

        %ADDING THE COMMENTED EQUATIONS SIGNIFICANTLY WORSENS THE
        %CONVERGENCE!!!!
        dE(3) = (E(4)-E(3))/dx;     %IF DON'T IMPOSE THIS BC, IT gives completely wrong results
        %dE(4) = (E(5)-E(4))/dx;  
        %dE(5) = (E(6)-E(5))/dx;  
        dE(nx-2) = (E(nx-2)-E(nx-3))/dx;
        %dE(nx-3) = (E(nx-3)-E(nx-4))/dx;
        %dE(nx-4) = (E(nx-4)-E(nx-5))/dx;
        %dE(nx-1) = 0;  %this isn't used anyway, so doesn't matter
       
        %Solve for new p
        old_p = p;
        for i = 3:nx-3        %only solve for the points inside the boundaries!          
            %attempt upwind standard derivatives for entire region
            %dE(i) = (E(i+1)-E(i))/dx;       
  
            p(i+1) = p(i) + dx*(-p(i)/E(i))*dE(i);   %there's divide by 0 issue here if E(i) = 0
   
            %stop run if NaN
            if isnan(p(i))
                stopstatement
            end
        end
        
        %weighting
        p = p*w + old_p*(1-w);
        
        error_p = max(abs(p-old_p)/abs(old_p))
    
       iter =  iter+1    
       
       if(Va == Va_max) %only for last run
          E_solution(iter) = E(46);    % save E at random point for each iter for convergence analysis
       end
    end
    
    %Calculate Jp
    for i = 3:nx-2 
        Jp(i) =  q*p_mob*p(i)*E(i);
    end
    
    %Save data
    str = sprintf('%.2f',Va);
    filename = [str 'V.txt'] 
    fid = fopen(fullfile('C:\Users\Tim\Documents\Duxbury group research\WENO\Mott-Gurney_law_WENO\Benchmarks\',filename),'w');   %w: Open or create new file for writing
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
V_final(Va_cnt) = 0;
    for i = 3:nx-2
        V_final(Va_cnt) = V_final(Va_cnt) + E(i)*dx;
    end  
    V_final

%% Final Plots

%Analytic Result Calculation
for i=3:nx-3
    E_theory1(i) = sqrt(2*x(i)*Jp(nx-3)/(epsilon*p_mob));
    E_theory2(i)= sqrt(2*x(i)*Jp(i)/(epsilon*p_mob));
    
%     E_theory1(i) = sqrt(2*x(i)*Jp(nx-3)/(epsilon*p_mob)+Ea^2);
%     E_theory2(i)= sqrt(2*x(i)*Jp(i)/(epsilon*p_mob)+Ea^2);
end

str = sprintf('%.2g', Va);

 h1 = plot(x(3:num_cell),p(3:num_cell));
 hold on
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);
 
 figure;
 h2 = plot(x(3:num_cell),E(3:num_cell));
 hold on
 plot(x(3:num_cell),E_theory1(3:num_cell));
 plot(x(3:num_cell),E_theory2(3:num_cell));
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);
 
 
 figure;
 h3 = plot(x(3:num_cell),Jp(3:num_cell));
 hold on
 title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
 
%  %JV curve
%  figure
%  h4 = plot(V,Jp_final);
%  hold on
%  plot(V, Jp_theory);
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
