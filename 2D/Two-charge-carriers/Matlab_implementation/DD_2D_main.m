%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         2D Drift Diffusion for Electrons and Holes with Finite Differences
%
%             -Scharfetter-Gummel discretization
%             -decoupled Gummel iterations method
%
%              Created by: Timofey Golubev (2018.05.29)
%
%
%     This includes the 2D poisson equation and 2D continuity/drift-diffusion
%     equations using Scharfetter-Gummel discretization. The Poisson equation
%     is solved first, and the solution of potential is used to calculate the
%     Bernoulli functions and solve the continuity eqn's.
%
%   Boundary conditions for Poisson equation are:
%
%     -a fixed voltage at (x,0) and (x, Nz) defined by V_bottomBC
%      and V_topBC which are defining the  electrodes
%    -insulating boundary conditions: V(0,z) = V(1,z) and
%     V(0,N+1) = V(1,N) (N is the last INTERIOR mesh point).
%     so the potential at the boundary is assumed to be the same as just inside
%     the boundary. Gradient of potential normal to these boundaries is 0.
%
%   Matrix equations are AV*V = bV, Ap*p = bp, and An*n = bn where AV, Ap, and An are sparse matrices
%   (generated using spdiag), for the Poisson and continuity equations.
%   V is the solution for electric potential, p is the solution for hole
%   density, n is solution for electron density
%   bV is the rhs of Poisson eqn which contains the charge densities and boundary conditions
%   bp is the rhs of hole continuity eqn which contains net generation rate
%   and BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

global num_cell N num_elements l_HTL_int Vt HTL_int_VBstep phi_a l_ETL_int ETL_int_VBstep phi_c BCP_int_VBstep N_dos p_topBC p_leftBC p_rightBC p_bottomBC Cn CV Cp n_leftBC  n_rightBC
global n_bottomBC n_topBC

%% Physical Constants
q =  1.60217646*10^-19;         %elementary charge, C
kb = 1.3806503*10^-23;          %Boltzmann const., J/k
T = 296.;                      %temperature:  Koster says they use room temperature...: I find that btw. 293 and 300 get no effect on JV curve...
epsilon_0 =  8.85418782*10^-12; %F/m
Vt = (kb*T)/q;

%% Simulation Setupt
Va_min = -0.5;               %volts
Va_max = 1.2;
increment = 0.01;         %for increasing V
num_V = floor((Va_max-Va_min)/increment)+1;   %number of V points

%Simulation parameters
w_eq = 0.01;               %For the 1st iteration (equilibrium run) to converge need a small weighting factor
w = 0.2;                 %set up of weighting factor
w_i = 0.2;
tolerance = 10^-12;        %error tolerance
tolerance_i =  5*10^-12;

%% System Setup
L = 10e-9;      %device thickness in meters
dx = 1e-9;                        %mesh size
num_cell = floor(L/dx);
N = num_cell -1;   %number of INTERIOR mesh points
num_elements = N^2;  %NOTE: this will specify number of elements in the solution vector V which = (num_cell +1 -2)^2 b/c we are in 2D
%(num_cell +1) = # of mesh pts, b/c matlab starts indexing from 1, then -2
%the endpts

%Electronic density of states of holes and electrons
N_VB = 10^24;         %density of states in valence band (holes)
N_CB = 10^24;         %density of states in conduction bands (electrons)
E_gap = 1.5;          %bandgap (in eV)
N_dos = 10^24.;            %scaling factor helps CV be on order of 1

%injection barriers
inj_a = 0.2;	%at anode
inj_c = 0.1;	%at cathode

%work functions of anode and cathode
WF_anode = 4.8;     
WF_cathode = 3.7;

Vbi = WF_anode - WF_cathode +inj_a +inj_c;  %built-in field


%% Define matrices of system parameters 

% Define relative dielectric constant matrix
%NOTE: epsilons will be defined at 1/2 integer points, so epsilons inside
%the cells, not at cell boundaries
%will use indexing: i + 1/2 is defined as i+1 for the index
epsilon = zeros(num_cell+1, num_cell +1);

%NOTE: I unfortunately can't define epsilon(0,..) in matlab, so the endpts
%are at 1...
for i = 1:num_cell+1
    epsilon(i,:) = 3.0;
end

% Define mobilities matrix
%can be position dependent
%using indexing: i+1/2 is defined as i+1 for index, just like the epsilons

for i = 1:num_cell+1
    for j = 1:num_cell+1
        p_mob(i,j) = 4.5*10^-6;
        n_mob(i,j) = 4.5*10^-6;
    end
end

mobil = 5.*10^-6;                    %scaling for mobility (for numerical stability when solving matrix equations)

p_mob = p_mob./mobil;
n_mob = n_mob./mobil;

%% Define Poisson equation boundary conditions and initial conditions
% Initial conditions
V_bottomBC = -((Vbi)/(2*Vt)-inj_a/Vt);
V_topBC = (Vbi)/(2*Vt)-inj_c/Vt;
diff = (V_topBC - V_bottomBC)/num_cell;
index = 0;
for j = 1:N  %corresponds to z coord
    index = index +1;
    V(index) = diff*j;
    for i = 2:N  %elements along the x direction assumed to have same V
        index = index +1;
        V(index) = V(index-1);
    end
end

%Side BCs will be filled in from V, since are insulating BC's
%i's in these BC's correspond to the x-value (z values are along a line,
%top and bottom)
%THESE NEED TO BE UPDATED AT EVERY ITERATION OF POISSON SOLVE
for i = 1:N
    V_leftBC(i) = V((i-1)*N + 1);  %just corresponds to values of V in the 1st subblock
    V_rightBC(i) = V(i*N);
end

%% Define continuity equn boundary and initial conditions
%these are scaled
n_bottomBC = N_CB*exp(-(E_gap-inj_a)/Vt)/N_dos;
p_bottomBC = N_VB*exp(-inj_a/Vt)/N_dos;
n_topBC = N_CB*exp(-inj_c/Vt)/N_dos;
p_topBC = N_VB*exp(-(E_gap-inj_c)/Vt)/N_dos;

%define initial conditions as min value of BCs
min_dense = min(n_bottomBC, p_topBC);
for i = 1:num_elements
    p(i,1) = min_dense;
    n(i,1) = min_dense;
end

for i = 1:N
    p_leftBC(i) = p((i-1)*N + 1);  %just corresponds to values in the 1st subblock
    p_rightBC(i) = p(i*N);
    n_leftBC(i) = n((i-1)*N + 1);  %just corresponds to values in the 1st subblock
    n_rightBC(i)= n(i*N);
end

%form matrices for easy filling of bV
n_matrix = reshape(n,N,N);
p_matrix = reshape(p,N,N);

% Set up Poisson matrix equation
AV = SetAV_2D(epsilon);
%spy(AV);  %allows to see matrix structure, very useful!

%% Solver Loop
%Scaling coefficients
Cp = dx^2/(Vt*N*mobil);          %note: scaled p_mob and n_mob are inside matrices
Cn = dx^2/(Vt*N*mobil);
CV = N*dx^2*q/(epsilon_0*Vt);    %relative permitivity was moved into the matrix

%for now, have no generation rate
Up = zeros(N,N);
Un = Up;

iter = 1;
error_np =  1;
while(error_np > tolerance)

    %set up rhs of Poisson equation. Note for epsilons, are assuming that
    %epsilons at the boundaries are the same as espilon 1 cell into interior of
    %device
    index = 0;
    for j = 1:N
        if(j ==1)  %different for 1st subblock
            for i = 1:N
                index = index +1;
                if (i==1)  %1st element has 2 BC's
                    bV(index,1) = CV*(p_matrix(i,j)-n_matrix(i,j)) + epsilon(i,j)*(V_leftBC(1) + V_bottomBC);   %RECALL matrix has +2 down diagonals, sign flipped from 1D version
                elseif (i==N)
                    bV(index,1) = CV*(p_matrix(i,j)-n_matrix(i,j)) + epsilon(i,j)*(V_rightBC(1) + V_bottomBC);
                else
                    bV(index,1) = CV*(p_matrix(i,j)-n_matrix(i,j)) + epsilon(i,j)*V_bottomBC;
                end
            end
        elseif(j == N)  %different for last subblock
            for i = 1:N
                index = index +1;
                if (i==1)  %1st element has 2 BC's
                    bV(index,1) = CV*(p_matrix(i,j)-n_matrix(i,j)) + epsilon(i,j)*(V_leftBC(N) + V_topBC);
                elseif (i==N)
                    bV(index,1) = CV*(p_matrix(i,j)-n_matrix(i,j)) + epsilon(i,j)*(V_rightBC(N) + V_topBC);
                else
                    bV(index,1) = CV*(p_matrix(i,j)-n_matrix(i,j)) + epsilon(i,j)*V_topBC;
                end
            end
        else %interior subblocks
            for i = 1:N
                index = index +1;
                if(i==1)
                    bV(index,1) = CV*(p_matrix(i,j)-n_matrix(i,j)) + epsilon(i,j)*V_leftBC(j);
                elseif(i==N)
                    bV(index,1) = CV*(p_matrix(i,j)-n_matrix(i,j)) + epsilon(i,j)*V_rightBC(j);
                else
                    bV(index,1) = CV*(p_matrix(i,j)-n_matrix(i,j));
                end
            end
        end
    end

    %solve for V
    oldV = V;
    newV = AV\bV;
    
    if(iter >1)
        V = newV*w + oldV*(1.-w);
    else
        V = newV;  %no mixing for 1st iter
    end
    
    %Update side boundary conditions
    %THESE NEED TO BE UPDATED AT EVERY ITERATION OF POISSON SOLVE
    for i = 1:N
        V_leftBC(i) = V((i-1)*N + 1);  %just corresponds to values of V in the 1st subblock
        V_rightBC(i) = V(i*N);
    end

    %reshape the V vector into the V matrix
    V_matrix = reshape(V,N,N);
    %add on the BC's to get full potential matrix
    fullV(2:N+1,2:N+1) = V_matrix;
    fullV(:,1) = V_bottomBC;
    fullV(:,N+2) = V_topBC;
    fullV(1,2:N+1) = V_leftBC;
    fullV(N+2,2:N+1) = V_rightBC;
    fullV(N+2, N+2) = V_topBC;  %assume that this far exterior corner has same V as rest of the top


    %% Set up continuity equation matrix
    Bernoulli_p_values = Calculate_Bernoullis_p(fullV); %the values are returned as a struct
    Bernoulli_n_values = Calculate_Bernoullis_n(fullV);
    Ap = SetAp_2D(p_mob, Bernoulli_p_values);           %send bernoulli values struct to matrix setup function
    An = SetAn_2D(n_mob, Bernoulli_n_values);
    bp = Setbp(Bernoulli_p_values, p_mob, Up);
    bn = Setbn(Bernoulli_n_values, n_mob, Un);
    
    oldp = p;
    newp = Ap\bp;
    
    oldn = n;
    newn = An\bn;
    
    old_error =  error_np;
    count = 0;
    error_np_matrix = zeros(1,num_elements); %need to reset the matrix b/c when more newp's become 0, need to remove that error element from matrix.
    for i = 1:num_elements  %recall that newp, oldp are stored as vectors
            if(newp(i) ~=0 && newn(i) ~=0)
                count = count+1;  %counts number of non zero error calculations
                error_np_matrix(count) = (abs(newp(i)-oldp(i)) + abs(newn(i) - oldn(i)))/abs(oldp(i) + oldn(i));  %need the dot slash, otherwise it tries to do matrix operation! %ERROR SHOULD BE CALCULATED BEFORE WEIGHTING
            end
    end
    error_np = max(error_np_matrix)
    
    w
    tolerance
    
    %weighting
    p = newp*w + oldp*(1.-w);
    n = newn*w + oldn*(1.-w);
    
    % update bc's and matrix of hole densities--------------------------------------------------------------
    %Update side boundary conditions
    %THESE NEED TO BE UPDATED AT EVERY ITERATION OF POISSON SOLVE
    for i = 1:N
        p_leftBC(i) = p((i-1)*N + 1);  %just corresponds to values of V in the 1st subblock
        p_rightBC(i)= p(i*N);
    end

    %reshape the p vector into the p matrix
    p_matrix = reshape(p,N,N);
    %add on the BC's to get full potential matrix
    fullp(2:N+1,2:N+1) = p_matrix;
    fullp(:,1) = p_bottomBC;
    fullp(:,N+2) = p_topBC;
    fullp(1,2:N+1) = p_leftBC;
    fullp(N+2,2:N+1) = p_rightBC;
    fullp(N+2, N+2) = p_topBC;  %assume that this far exterior corner has same V as rest of the top

    % update bc's and matrix of electron densities ---------------------------------------------------------
    
     for i = 1:N
        n_leftBC(i) = n((i-1)*N + 1);  %just corresponds to values of V in the 1st subblock
        n_rightBC(i)= n(i*N);
    end

    %reshape then vector into the n matrix
    n_matrix = reshape(n,N,N);
    %add on the BC's to get full potential matrix
    fulln(2:N+1,2:N+1) = n_matrix;
    fulln(:,1) = n_bottomBC;
    fulln(:,N+2) = n_topBC;
    fulln(1,2:N+1) =n_leftBC;
    fulln(N+2,2:N+1) = n_rightBC;
    fulln(N+2, N+2) = n_topBC;  %assume that this far exterior corner has same V as rest of the top
    
    iter = iter +1

end

%plot V
surf(1:N+2,1:N+2, Vt*fullV)

figure
%plot p
surf(1:N+2,1:N+2, log(N*fullp))  

hold on
surf(1:N+2,1:N+2, log(N*fulln))  



