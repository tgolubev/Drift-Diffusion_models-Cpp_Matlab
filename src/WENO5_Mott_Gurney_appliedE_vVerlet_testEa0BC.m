%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving Poisson + Mott-Gurney limit system (for holes) with 5th order
%             Weighted Essentially Non-Oscilatory (WENO5)
%
%                        dE/dx = (q/epsilon)*p
%                E*dp/dx + p*dE/dx = 0  for x in [0,L]
%                  
%          This version directly applies an E to one electrode.
%       Uses modified velocity Verlet (with WENO) for Poisson eqn.
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
L = 100*10^-9;             %device length in meters
num_cell = 5000;            % number of cells
p_initial =  10^27;        %initial hole density
p_mob = 2.0*10^-8;         %hole mobility

Va_min = 190;              %volts
Va_max = 200;    
V_increment = 1;           %for increasing V
Ea_min = Va_min/L;         %V/m
Ea_max = Va_max/L;         %maximum applied E
%increment = 0.1*10^8;%V_increment/L; %for increasing E

%Simulation parameters
w= 1.;      %weighting factor
constant_p_i = false;
tolerance = 10^-14;   %error tolerance       
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
%Note: p(1), p(2),p(nx-1), and p(nx) values are not needed b/c they are never used (by
%default are 0).

if(constant_p_i)
    for i = 3:nx-2
        p(i) = p_initial;
    end
else
    %linearly decreasing p: THIS DOESN'T WORK as well
    p(3) = p_initial;
    for i = 3:nx-3     
        dp = p_initial/(num_cell+1);
        p(i+1) = p(i)-dp;
    end
end

%initial condition: make E be constant
for i = 3:nx-2
    E(i) = Ea_min;
end

Ea_cnt = 0;    
for Va = Va_min:V_increment:Va_max
    
    Ea_cnt = Ea_cnt+1;
    
    %Boundary conditions
    %because need 2 ghost points to left and right and matlab indices arrays from 1, means that E(x=0) =
    %E(i=3). And E(nx-2) = E(x=L)
    %NOTE: IF SET THESE TO = Ea , then have convergence issues
    E(1) = 0;
    E(2) = 0;
    %E(nx-1) = 0;
    %E(nx) = 0;
    
    
    %CAN'T FORCE E(3) to be 0: b/c will result in blowup issues later: so
    %need to set E(3) to equal (V(4)-V(3))/dx: see below
    
    E(3) = 0;   %ONLY INITIALLY WE ASSUME THIS, THEN WE RESET E
    
    %E(1) = E(3);
    %E(2) = E(3);
    E(nx) = E(nx-2);
    E(nx-1) = E(nx-2);
    
  
   
    %these are needed for using WENO to find dp/dx
    %NOTE: IF SET THESE TO = P_INITIAL , then have convergence issues
    p(1) = 0;
    p(2) = 0;
    p(nx) = 0;
    p(nx-1) = 0;
    
    %right side boundary will be determined by Newton's method solving of
    %Poisson eqn.

    %timing
    tic
    
    %% Solver Loop
    iter = 1;
    error_p =  1.0;
    while error_p > tolerance
        
        %Poisson equation with velocity Verlet method
        
        %use WENO to estimate dp/dx
        dp = residual(p,flux,dflux,dx,nx,fluxsplit);
        
        dp(3) = (p(4)-p(3))/dx;
        dp(nx-2) = (p(nx-2)-p(nx-3))/dx;
       
        p_calc(3) = p(2) + dx*dp(3);    %0.5*dx*(dp(3) + dp(2)); if use WENO this averaging will be wrong b/c dp(2) = 0, so taking only 1/2 of the slope
        
        for i = 3:nx-3       %only solve for E inside the device (points 1,2, nx-1, and nx are outside device. 
            %dp(i) = (p(i+1)-p(i))/dx;  %try using finite differences
            %E(i+1) = E(i) + (q/epsilon)*p(i);    %Euler method
            p_calc(i+1) = p(i) + dx*dp(i); %  0.5*dx*(dp(i+1) + dp(i));  %not necessary to do the 1/2 averaging if use WENO since our dp(i) is calculated more sophisticated.
            E(i+1) = E(i) + (q/epsilon)*p_calc(i)*dx+ 0.5*dx^2*(q/epsilon)*dp(i);   
        end
        
        % need to set E(3) by using that V = -integral (E*dx)
        %so V(4) = integragral (E*dx) from i = 3 to 4
        %and V(3) we know is Va
        % then use E(3) = (V(4)-V(3))/dx
        %assume initially  E(3) = 0 as set above
        E(3) = -(E(4)*dx-Va)/dx;
     sdfsdf   
        
       %integrate E to enforce the applied V condition
        V_actual = 0;
        for i = 3:nx-2
            V_actual = V_actual + E(i)*dx;
        end 
        E = E*Va/V_actual  %normalizes E so have correct applied voltage
    
        
        %weight E
        %old_E = E;
        %E = w*E + (1.-w)*old_E;
        
        %dE = weno approx for dE/dx
        dE = residual(E,flux,dflux,dx,nx,fluxsplit);     %this calculates the entire dE array
        
        %RIGHT NOW HAVING ISSUE THAT THE dE(3) is too large: having
        %issue doing the boundary derivative estimate properly
        
        
        %adding the BELOW COMMENTED EQUATIONS ACTUALLY MESSES UP THE
        %CONVERGENCE!!
        dE(3) = (E(4)-E(3))/dx;     %IF DON'T IMPOSE THIS BC, IT REALLY FAILS!
        %dE(4) = (E(5)-E(4))/dx; 
        %dE(5) = (E(6)-E(5))/dx;  
        dE(nx-2) = (E(nx-2)-E(nx-3))/dx;
        %dE(nx-3) = (E(nx-3)-E(nx-4))/dx;
      
        % Solve for new p
        old_p = p;    %for computing error
        for i = 3:nx-4        %only solve for the points inside the boundaries!  
          
            %attempt upwind standard derivatives for entire region
            %dE(i) = (E(i+1)-E(i))/dx;       %THIS WORKS but takes more iterations 
            %for convergence and approach to convergence is oscillatory
  
            p(i+1) = p(i) + dx*(-p(i)/E(i))*dE(i);   %there's divide by 0 issue here if E(i) = 0
             
            %stop run if NaN
            if isnan(p(i))
                stopstatement
            end
        end
        
        %try weighting
        p = w*p + (1.-w)*old_p;
        
        error_p = max(abs(p-old_p)/abs(old_p))
        
        %adjust BC's for p: IS NOT NECESSARY B/C p(3) is never recalculated
        %p(3) = p_initial;   %for conservation of particles
    
       iter =  iter+1   
       
       if(Va == Va_max) %only for last run
            E_solution(iter) = E(50);    % save E at random point for each iter for convergence analysis
       end
    end
    
    %Calculate Jp
    for i = 3:nx-2 
        Jp(i) =  q*p_mob*p(i)*E(i);
    end
  
    %Ea
    
    %Calculations for JV curve
    %Calculate Va_final for each Ea by integrating E
    V(Ea_cnt) = 0;
    Jp_final(Ea_cnt) = Jp(nx-3);  %just pick Jp at the right side
    for i = 3:nx-2
        V(Ea_cnt) = V(Ea_cnt) + E(i)*dx;
    end    
    
    %sqrt(2Jx/p_mob*epsilon)). THIS GIVES HUGE NUMBERS!
    %Va_final(Ea_cnt) = (2/3)*sqrt(2*Jp(nx-3)/(epsilon*p_mob))*L^(2/3)
    
    % Save data
    %str = sprintf('%.2f',Ea*L);
    str = sprintf('%.2g',Va);
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


%% Final Plots

%Analytic Result Calculation
for i=3:nx-3
    E_theory1(i) = sqrt(2*x(i)*Jp(nx-3)/(epsilon*p_mob));
    E_theory2(i)= sqrt(2*x(i)*Jp(i)/(epsilon*p_mob));
%     E_theory1(i) = sqrt(2*x(i)*Jp(nx-3)/(epsilon*p_mob)+Ea^2);
%     E_theory2(i)= sqrt(2*x(i)*Jp(i)/(epsilon*p_mob)+Ea^2);
end
for j = 1:Ea_cnt
    Jp_theory(j) = (9*p_mob*epsilon*(V(j))^2)/(8*L^3);
end

%Va_final = Ea*L;
%str=sprintf('%.3f', Va_final);   %.3 precision operator sets 3 decimals
str = sprintf('%.2g', Ea);

 h1 = plot(x(3:num_cell),p(3:num_cell));
 hold on
 title(['Ea =', str, 'V/m'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);
 
 figure;
 h2 = plot(x(3:num_cell),E(3:num_cell));
 hold on
 plot(x(3:num_cell),E_theory1(3:num_cell));
 plot(x(3:num_cell),E_theory2(3:num_cell));
 title(['Ea =', str, 'V/m'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);
 
 
 figure;
 h3 = plot(x(3:num_cell),Jp(3:num_cell));
 hold on
 title(['Ea =', str, 'V/m'],'interpreter','latex','FontSize',16);
 xlabel('Position ($m$)','interpreter','latex','FontSize',14);
 ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
 
 %JV curve
 figure
 h4 = plot(V,Jp_final);
 hold on
 plot(V, Jp_theory);
 xlabel('Voltage (V)','interpreter','latex','FontSize',14);
 ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);
 
 %convergence analysis
 iterations = 1:iter;
 figure
 h5 = plot(iterations, E_solution);
 title(['E convergence', str, 'V/m'],'interpreter','latex','FontSize',16);
 xlabel('Iterations','interpreter','latex','FontSize',14);
 ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);
 
 hold off

