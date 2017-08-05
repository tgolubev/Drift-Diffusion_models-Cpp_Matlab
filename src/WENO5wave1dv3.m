%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Solving Poisson + Mott-Gurney limit system with 5th order
%             Weighted Essentially Non-Oscilaroty (WENO5)
%
%                E*dp/dx + p*dE/dx = 0  for x \in [0,L]
%                  where f = f(u): linear/nonlinear
%
%            
%       I modified this code based on original code of 1D wave eqn
%                 by Manuel Diaz, manuel.ade'at'gmail.com 
%              Institute of Applied Mechanics, 2012.08.20
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: Jiang & Shu; Efficient Implementation of Weighted ENO Schemes
% JCP. vol 126, 202-228 (1996)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Finite Difference Implementation of WENO5 with SSP-RK33 time
% integration method. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
      cfl = 0.30;   % Courant Number
      L = 70;      %device length 
      nx = 150;    % number of cells
     tEnd = 2*pi;   % End time
  BC_type = 3;      % {1} Dirichlet, {2} Neumann, {3} Periodic
fluxsplit = 2;      % {1} Godunov, {2} Global LF, {3} Local LF

%% Define our Flux function
fluxtype = 'linear';

switch fluxtype
    case 'linear'
        c=1; flux = @(w) c*w;
        dflux = @(w) c*ones(size(w));
    case 'nonlinear' % Burgers'
        flux = @(w) w.^2/2;
        dflux = @(w) w;
end

% Domain Discretization
a=0; b=L; x0=linspace(a,b,nx); dx=(b-a)/nx;   %x0 is positions array
% Depending on the degree of the WENO stencil we are using, 2 more ghost
% will have to be added to the final domain.
%x is the array of 5 points which will be used to estimate dJ/dx
x= [-2*dx,-dx,x0,b+dx,b+2*dx]; nx = length(x);  

% Initial Condition
u0 = IC(x,6);

% Plot range
domain = 3:nx-2; % don't show gost cells
plotrange = [a,b,min(u0)-0.1,max(u0)+0.1];

%% Solver Loop
% Beause the max slope, f'(u) = u, may change as the time steps progress
% we cannot fix the size of the time steps. Therefore we need to recompute
% the size the next time step, dt, at the beginning of each iteration to
% ensure stability.

% load initial conditions
t=0; it=0; u=u0;

while t < tEnd
    uo = u;
    
    % Update time
    dt = cfl*dx/abs(max(u(domain))); t=t+dt;

    % iteration counter
    it=it+1;  
    
    % 1st stage: dF = weno approx for df/dx
    dF = residual(u,flux,dflux,dx,nx,fluxsplit);
    u = uo-dt*dF;
    
    % 2nd Stage
    dF = residual(u,flux,dflux,dx,nx,fluxsplit);
    u = 0.75*uo+0.25*(u-dt*dF);

    % 3rd stage
    dF = residual(u,flux,dflux,dx,nx,fluxsplit);
    u = (uo+2*(u-dt*dF))/3;
    
    % Plot solution
    plot(x,u0,'-x',x(domain),u(domain),'.'); axis(plotrange)
    
    %if rem(it,10) == 0
        drawnow;
    %end
end
%% Final Plot
plot(x,u0,'-x',x(domain),u(domain),'-'); 
axis([a,b,min(u0)-0.1,max(u0)+0.1])
title('WENO5, Cell Averages plot','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);