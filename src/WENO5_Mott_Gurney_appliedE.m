%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving Poisson + Mott-Gurney limit system (for holes) with 5th order
%             Weighted Essentially Non-Oscilatory (WENO5)
%
%                E*dp/dx + p*dE/dx = 0  for x in [0,L]
%                  
%          This version directly applies an E to one electrode.
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
num_cell = 500;            % number of cells
p_initial =  10^27;        %initial hole density
p_mob = 2.0*10^-8;         %hole mobility

Va_min = 190;              %volts
Va_max = 200;    
V_increment = 1;           %for increasing V
Ea_min = Va_min/L;         %V/m
Ea_max = Va_max/L;         %maximum applied E
increment = V_increment/L; %for increasing E

%Simulation parameters
constant_p_i = true;
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

%Doesn't seem to make difference whether use constant_p_i or linearlrly
%decreasing
if(constant_p_i)
    for i = 3:nx-2
        p(i) = p_initial;
    end
else
    %linearly decreasing p
    p(3) = p_initial;
    for i = 3:nx-3     
        dp = p_initial/(num_cell+1);
        p(i+1) = p(i)-dp;
    end
end

    
for Ea = Ea_min:increment:Ea_max
    
    %Boundary conditions
    %because need 2 ghost points to left and right and matlab indices arrays from 1, means that E(x=0) =
    %E(i=3). And E(nx-2) = E(x=L)
    E(1) = 0;
    E(2) = 0;
    E(3) = Ea;
    
    %right side boundary will be determined by Newton's method solving of
    %Poisson eqn.

    %% Solver Loop
    iter = 0;
    error_p =  1.0;
    while error_p > tolerance
        
        %Poisson equation with Newton's method
        for i = 3:nx-3       %only solve for E inside the device (points 1,2, nx-1, and nx are outside device. point 3 is used to enforce E(x=0) = 0 BC.
            E(i+1) = E(i) + (q/epsilon)*p(i)*dx;   %means with initial constant p(i), this will be linear
        end
        
        %Define E at right side ghost points
        E(nx-1) = E(nx-2);
        E(nx) = E(nx-2);
        
        %dE = weno approx for dE/dx
        dE = residual(E,flux,dflux,dx,nx,fluxsplit);     %this calculates the entire dE array
        
        %RIGHT NOW HAVING ISSUE THAT THE dE(3) is too large: having
        %issue doing the boundary derivative estimate properly

        %try for i=3 to estimate derivative by  standard finite
        %difference:  THIS SOLVED BOTH THE KINK IN E GRAPH AND ISSUE OF STARTING AT I=3 NOT WORKING!

        dE(3) = (E(4)-E(3))/dx;
        
        %Solve for new p
        for i = 3:nx-3        %only solve for the points inside the boundaries!  
          
            %attempt upwind standard derivatives for entire region
            %dE(i) = (E(i+1)-E(i))/dx;       %THIS WORKS!
           
            old_p = p;
  
            p(i+1) = p(i) + dx*(-p(i)/E(i))*dE(i);   %there's divide by 0 issue here if E(i) = 0
            
            error_p = max(abs(p-old_p)/abs(old_p));
         

            %stop run if NaN
            if isnan(p(i))
                stopstatement
            end

        end
        
        %adjust BC's for p
        p(1) = 0;
        p(2) = 0;
        p(3) = p_initial;   %for conservation of particles
        p(nx) = 0;
        p(nx-1) = 0;
    
       iter =  iter+1;    
       iter
    end
    
    %Calculate Jp
    for i = 3:nx-2 
        Jp(i) =  q*p_mob*p(i)*E(i);
    end
  
    %Ea
    
    %Save data
    str = sprintf('%.2f',Ea*L);
    filename = [str 'V.txt'] 
    fid = fopen(fullfile('C:\Users\Tim\Documents\Duxbury group research\WENO\Mott-Gurney_law_WENO\Benchmarks\',filename),'w');   %w: Open or create new file for writing
    %fullfile allows to make filename from parts
    for i = 3:nx-2
        fprintf(fid,'%.8e %.8e %.8e %.8e\r\n ',x(i), p(i), E(i), Jp(i)); 
        %f means scientific notation, \r\n: start new line for each value
        %for each new desired column, put its own % sign
      
    end
    fclose(fid);
    
end


%% Final Plots

%Analytic Result Calculation
for i=3:nx-3
    E_theory1(i) = sqrt(2*x(i)*Jp(nx-3)/(epsilon*p_mob));
    E_theory2(i)= sqrt(2*x(i)*Jp(i)/(epsilon*p_mob));
    
%     E_theory1(i) = sqrt(2*x(i)*Jp(nx-3)/(epsilon*p_mob)+Ea^2);
%     E_theory2(i)= sqrt(2*x(i)*Jp(i)/(epsilon*p_mob)+Ea^2);
end

Va_final = Ea*L;
str=sprintf('%.3f', Va_final);   %.3 precision operator sets 3 decimals

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

 
hold off

