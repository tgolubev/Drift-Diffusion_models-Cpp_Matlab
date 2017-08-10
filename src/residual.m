function dF = residual(u,flux,dflux,dx,nx,fluxsplit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Compute the Residual using WENO5
%
%                         Residual = dF/dx
%             where  F= is our WENO Flux across cells
%
%   This is the WENO finite diff. approximation of the 1st derivative
%
%  Slightly modified from original code by Manuel Diaz, NTU, 2012.08.20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flux Spliting
[vp,vn] = fluxSplit(u,flux,dflux,fluxsplit);

% Allocate arrays
hn = zeros(size(vn)); hp = zeros(size(vp));


% Reconstruc h+ and h- fluxes
for i = 3:nx-2      % from cell 3 to nx-2
    sr = i-2:i+2;   % Stencil range
    [hn(i),hp(i-1)] = WENO5_1d_reconstruction(vp(sr),vn(sr));    
end
% Total flux
h = hn + hp;

% Set BCs   
%h(1:2) = 0;     %to enforce conservation, h(2) is hi-1/2 for the x = 0 point. 
%h(nx-1:nx) = 0;   % is hi+1/2 for the x=L point: note the i+1/2 are defined as hn(i) and i-1/2 as hp(i-1)

% Set Periodic BC
%h(1:2) = h(nx-4:nx-3);
%h(nx-1:nx) = h(3:4);


% Formulate Left and Right fluxes, equiv: % h(i)-h(i-1)
h_right = h; h_left = circshift(h,[0,1]);

% residual
dF = (h_right - h_left)/dx;

