%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         3D Drift Diffusion for Holes with Finite Differences

%           THIS IS THE SINGLE CARRIER VERSION--> NO GENERATION 
%       UNIFORM TOP BC
%
%             -Scharfetter-Gummel discretization
%             -decoupled Gummel iterations method
%
%            Created by: Timofey Golubev (2018.07.06)
%
%     This includes the 3D poisson equation and 3D continuity/drift-diffusion
%     equations using Scharfetter-Gummel discretization. The Poisson equation
%     is solved first, and the solution of potential is used to calculate the
%     Bernoulli functions and solve the continuity eqn's.
%
%   Boundary conditions for Poisson equation are:
%
%     -a fixed voltage at (x,0) and (x, Nz) defined by V_bottomBC
%      and V_topBC which are defining the  electrodes
%
%    -insulating boundary conditions: V(0,y,z) = V(1,y,z) and
%     V(N+1,y,z) = V(N,y,z) (N is the last INTERIOR mesh point).
%     so the potential at the boundary is assumed to be the same as just inside
%     the boundary. Gradient of potential normal to these boundaries is 0.
%    V(x,0,z) = V(x,1,z) and V(x,N+1,z) = V(x,N,z)
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
global num_cell_x num_cell_y num_cell_z Nx Ny Nz  num_elements Vt N_dos p_topBC p_bottomBC Cn CV Cp 
global V_bottomBC V_topBC  dx dy  dz

%% Physical Constants

q =  1.60217646*10^-19;         %elementary charge (C)
kb = 1.3806503*10^-23;          %Boltzmann const. (J/k)
T = 296.;                       %temperature (K)
epsilon_0 =  8.85418782*10^-12; %F/m
Vt = (kb*T)/q;

%% Simulation Setup

%Voltage sweep loop
Va_min = 0.1;            %volts
Va_max = 10.0;
increment = 0.1;         %by which to increase V
num_V = floor((Va_max-Va_min)/increment)+1;   %number of V points

%Simulation parameters
w_eq = 0.2;               %linear mixing factor for 1st convergence (0 applied voltage, no generation equilibrium case)
w_i = 0.2;                 %starting linear mixing factor for Va_min (will be auto decreased if convergence not reached)
tolerance = 5*10^-12;        %error tolerance
tolerance_i =  5*10^-12;     %initial error tolerance, will be increased if can't converge to this level

%% System Setup

Lx = 6.0000001e-9;     %there's some integer rounding issue, so use this .0000001
Ly = 6.0000001e-9;
Lz = 6.0000001e-9;   
dx = 2e-9;                        %mesh size
dy = 2e-9;  
dz = 2e-9;  
num_cell_x = floor(Lx/dx);
num_cell_y = floor(Ly/dy);
num_cell_z = floor(Lz/dz);
Nx = num_cell_x -1;       %number of INTERIOR mesh points (total mesh pts = num_cell +1 b/c matlab indixes from 1)
Ny = num_cell_y -1; 
Nz = num_cell_z -1; 
num_elements = (Nx+1)*(Ny+1)*(Nz+1);  %NOTE: this will specify number of elements in the solution vector V

% %Electronic density of states of holes 
 N_dos = 1.25*10^27.;       %scaling factor helps CV be on order of 1


%% Define matrices of system parameters

startProgram = tic;

%Preallocate vectors and matrices
fullV = zeros(Nx+2, Ny+2, Nz+2);
fullp = zeros(Nx+2, Ny+2, Nz+2);
% fulln = zeros(N+2, N+2, N+2);
Jp_Z = zeros(num_cell_x, num_cell_y, num_cell_z);
Jn_Z = zeros(num_cell_x, num_cell_y, num_cell_z);
Jp_X = zeros(num_cell_x, num_cell_y, num_cell_z);
% Jn_X = zeros(num_cell, num_cell, num_cell);
% Jp_Y = zeros(num_cell, num_cell, num_cell);
% Jn_Y = zeros(num_cell, num_cell, num_cell);
V_values = zeros(num_V+1,1);
J_total_Z_middle = zeros(num_V+1,1);

% Relative dielectric constant matrix (can be position dependent)
%Epsilons are defined at 1/2 integer points, so epsilons inside
%the cells, not at cell boundaries
%will use indexing: i + 1/2 is defined as i+1 for the index
epsilon = 3.0*ones(num_cell_x+2, num_cell_y+2, num_cell_z+2);
for i =  1:num_cell_x
    for j = 1:num_cell_y
        for k = 1:num_cell_z
            epsilon_avged.eps_X_avg(i+1,j+1,k+1) = (epsilon(i,j,k) + epsilon(i,j+1,k) + epsilon(i,j,k+1) + epsilon(i,j+1,k+1))./4.;
            epsilon_avged.eps_Y_avg(i+1,j+1,k+1) = (epsilon(i,j,k) + epsilon(i+1,j,k) + epsilon(i,j,k+1) + epsilon(i+1,j,k+1))./4.;
            epsilon_avged.eps_Z_avg(i+1,j+1,k+1) = (epsilon(i,j,k) + epsilon(i+1,j,k) + epsilon(i,j+1,k) + epsilon(i+1,j+1,k))./4.;

            %add num_cell+2 values for i to account for extra bndry pt
            epsilon_avged.eps_X_avg(num_cell_x+2,j+1,k+1) = epsilon_avged.eps_X_avg(2,j+1,k+1);  %2 b/c the epsilon avg indices start from 2 (like Bernoullis)
            epsilon_avged.eps_Y_avg(num_cell_x+2,j+1,k+1) = epsilon_avged.eps_Y_avg(2,j+1,k+1);
            epsilon_avged.eps_Z_avg(num_cell_x+2,j+1,k+1) = epsilon_avged.eps_Z_avg(2,j+1,k+1);
            
            %add num_cell+2 values for j to account for extra bndry pt
            epsilon_avged.eps_X_avg(i+1,num_cell_y+2,k+1) = epsilon_avged.eps_X_avg(i+1,2,k+1);
            epsilon_avged.eps_Y_avg(i+1,num_cell_y+2,k+1) = epsilon_avged.eps_Y_avg(i+1,2,k+1);
            epsilon_avged.eps_Z_avg(i+1,num_cell_y+2,k+1) = epsilon_avged.eps_Z_avg(i+1,2,k+1);
        end
       %NOTE: use inside the device value for the num_cell+2 values b/c
       %these actually correspond to num_cell values, since with  Neuman
       %BC's, we use the value inside the device twice, when discretizing
       %at the boundary.
        epsilon_avged.eps_X_avg(i+1,j+1,num_cell_z+2) = epsilon(i+1,j+1,k+1);        %assume the average epsilon at the bndry pt (top electrode) = same as epsilon just inside.
        epsilon_avged.eps_Y_avg(i+1,j+1,num_cell_z+2) = epsilon(i+1,j+1,k+1);
        epsilon_avged.eps_Z_avg(i+1,j+1,num_cell_z+2) = epsilon(i+1,j+1,k+1);
    end
end

%-----------------Mobilities setup-----------------------------------------
% Define mobilities matrix (can be position dependent)
p_mob = (10^-8)*ones(num_cell_x+2, num_cell_y+2, num_cell_z+2);
% n_mob = p_mob;

mobil = 10^-8;        %scaling for mobility 

p_mob = p_mob./mobil;
% n_mob = n_mob./mobil;

%Pre-calculate the mobility averages 
%using indexing: i+1/2 is defined as i+1 for index, just like the epsilons
p_mob_avged.p_mob_X_avg = zeros(num_cell_x+2, num_cell_y+2, num_cell_z+1);
p_mob_avged.p_mob_Y_avg = zeros(num_cell_x+2, num_cell_y+2, num_cell_z+1);
p_mob_avged.p_mob_Z_avg = zeros(num_cell_x+2, num_cell_y+2, num_cell_z+1);

% n_mob_avged.n_mob_X_avg = zeros(num_cell+1, num_cell+1, num_cell+1);
% n_mob_avged.n_mob_Y_avg = zeros(num_cell+1, num_cell+1, num_cell+1);
% n_mob_avged.n_mob_Z_avg = zeros(num_cell+1, num_cell+1, num_cell+1);

%NOTE: THIS SHOULD BE LABELED AS +1, so i.e. averaging pmob (1) with
%pmob(2) should be called pmobavged (2) -> then will correspond with how I
%define  the bernoulli fnc's!!

for i =  1:num_cell_x
    for j = 1:num_cell_y  %to account for extra bndry pt added at right side
        for k = 1:num_cell_z
            p_mob_avged.p_mob_X_avg(i+1,j+1,k+1) = (p_mob(i,j,k) + p_mob(i,j+1,k) + p_mob(i,j,k+1) + p_mob(i,j+1,k+1))./4.;
            p_mob_avged.p_mob_Y_avg(i+1,j+1,k+1) = (p_mob(i,j,k) + p_mob(i+1,j,k) + p_mob(i,j,k+1) + p_mob(i+1,j,k+1))./4.;
            p_mob_avged.p_mob_Z_avg(i+1,j+1,k+1) = (p_mob(i,j,k) + p_mob(i+1,j,k) + p_mob(i,j+1,k) + p_mob(i+1,j+1,k))./4.;
            
            %             n_mob_avged.n_mob_X_avg(i,j,k) = (n_mob(i,j,k) + n_mob(i,j+1,k) + n_mob(i,j,k+1) + n_mob(i,j+1,k+1))./4.;
            %             n_mob_avged.n_mob_Y_avg(i,j,k) = (n_mob(i,j,k) + n_mob(i+1,j,k) + n_mob(i,j,k+1) + n_mob(i+1,j,k+1))./4.;
            %             n_mob_avged.n_mob_Z_avg(i,j,k) = (n_mob(i,j,k) + n_mob(i+1,j,k) + n_mob(i,j+1,k) + n_mob(i+1,j+1,k))./4.;
            %add num_cell+2 values for i to account for extra bndry pt
            p_mob_avged.p_mob_X_avg(num_cell_x+2,j+1,k+1) = p_mob_avged.p_mob_X_avg(2,j+1,k+1);
            p_mob_avged.p_mob_Y_avg(num_cell_x+2,j+1,k+1) = p_mob_avged.p_mob_Y_avg(2,j+1,k+1);
            p_mob_avged.p_mob_Z_avg(num_cell_x+2,j+1,k+1) = p_mob_avged.p_mob_Z_avg(2,j+1,k+1);
            
            %add num_cell+2 values for j to account for extra bndry pt
            p_mob_avged.p_mob_X_avg(i+1,num_cell_y+2,k+1) = p_mob_avged.p_mob_X_avg(i+1,2,k+1);
            p_mob_avged.p_mob_Y_avg(i+1,num_cell_y+2,k+1) = p_mob_avged.p_mob_Y_avg(i+1,2,k+1);
            p_mob_avged.p_mob_Z_avg(i+1,num_cell_y+2,k+1) = p_mob_avged.p_mob_Z_avg(i+1,2,k+1);
            
        end
            %NOTE: use inside the device value for the num_cell+2 values b/c
       %these actually correspond to num_cell values, since with  Neuman
       %BC's, we use the value inside the device twice, when discretizing
       %at the boundary.
        p_mob_avged.p_mob_X_avg(i+1,j+1,num_cell_z+2) = p_mob_avged.p_mob_X_avg(i+1,j+1,k+1);
        p_mob_avged.p_mob_Y_avg(i+1,j+1,num_cell_z+2) = p_mob_avged.p_mob_Y_avg(i+1,j+1,k+1);
        p_mob_avged.p_mob_Z_avg(i+1,j+1,num_cell_z+2) = p_mob_avged.p_mob_Z_avg(i+1,j+1,k+1);
    end
end


%--------------------------------------------------------------------------
%Scaling coefficients
Cp = dz^2/(Vt*N_dos*mobil);          %note: scaled p_mob and n_mob are inside matrices
% Cn = dz^2/(Vt*N_dos*mobil);
CV = (N_dos*dz^2*q)/(epsilon_0*Vt);     %relative permitivity was moved into the matrix

%% Define Poisson equation boundary conditions and initial conditions

% Initial conditions
V_bottomBC(1:Nx+2, 1:Ny+2) = 0;  %needs to be matrix, so can add i.e. afm tip
V_topBC(1:Nx+2, 1:Ny+2) = Va_min/Vt;  %set to 0 here-->  will set the tip voltage = Va/Vt inside the Va loop
diff = (V_topBC(1,1) - V_bottomBC(1,1))/num_cell_z;
% V(1:N) = V_bottomBC + diff;  %define V's corresponding to 1st subblock here (1st interior row of system)
% index = 0;
V_matrix = zeros(Nx+1,Ny+1,Nz+1);

%SET UP V_matrix first--> then use permute and (:) to get V--> that will be easier
for k = 1:Nz+1
    for i = 1:Nx+1
        for j = 1:Ny+1
            V_matrix(i, j, k) = V_bottomBC(i+1,j+1) +  diff*k;  
        end
    end
end
permuted_V_matrix = permute(V_matrix, [3 2 1]);
V = permuted_V_matrix(:);


%-------------------------------------------------------------------------------------------------
%% Define continuity equn boundary and initial conditions
%these are scaled
% n_bottomBC = N_CB*exp(-(E_gap-inj_a)/Vt)/N_dos;
p_bottomBC(1:Nx+2,1:Ny+2) = 0;
% n_topBC = N_CB*exp(-inj_c/Vt)/N_dos;

%FOR DIRICHLET BC'S
p_topBC(1:Nx+2,1:Ny+2) = 1;


%THIS DOESN'T WORK!
% diff_p = (p_topBC(1,1) - p_bottomBC(1,1))/num_cell_z;
% p_guess_matrix = zeros(Nx+1,Ny+1,Nz+1);
% 

% %setup initial guess for p for continuity eqn solve
% for k = 1:Nz+1
%     for i = 1:Nx+1
%         for j = 1:Ny+1
%             p_guess_matrix(i, j, k) = p_bottomBC(i+1,j+1) +  diff_p*k; 
%         end
%     end
% end
% permuted_p_matrix = permute(p_guess_matrix, [3 2 1]);
% p_guess = permuted_p_matrix(:);


%define initial conditions as min value of BCs --> THIS are initial
%conditions for Poisson eqn
min_dense = min(p_bottomBC(1,1), p_topBC(1,1));
p = min_dense*ones(num_elements, 1);
% p = n;

%form matrices for easy filling of bp
% n_matrix = reshape(n,N,N,N);
p_matrix = reshape(p,Nx+1,Ny+1,Nz+1); %this is a reshape before solving, so keep as x,y,z format


%-------------------------------------------------------------------------------------------------

% Set up Poisson matrix equation
AV = SetAV_3D(epsilon_avged);
% [L,U] = lu(AV);  %do and LU factorization here--> since Poisson matrix doesn't change
%this will significantly speed up backslash, on LU factorized matrix
%spy(AV);  %allows to see matrix structure, very useful!

%% Main voltage loop
Va_cnt = 1;
for Va_cnt = 1:num_V 
    startLoop = tic;
    not_converged = false;
    not_cnv_cnt = 0;
    
    %stop the calculation if tolerance becomes too high
    if(tolerance >10^-5)
        break
    end

    if(Va_cnt > 0)
        tolerance = tolerance_i;       %reset tolerance back
        w=w_i;  %I'm refreshing w every time...
        Up = zeros(num_elements,1);  %needed for bV setting
        %         G = GenerationRate();  %only use it once, since stays constant
        Va = Va_min+increment*(Va_cnt-1);     %set Va value
        
        %Va = Va_max-increment*(Va_cnt-1);    %decrease Va by increment in each iteration
    end
    
    %Voltage boundary conditions
    V_bottomBC(1:Nx+2, 1:Ny+2) = 0;
    %note: N+1 = num_cell --> these indices for tip positionare consistent with ones used
    %in Ap and AV set
    
    %for Dirichlet BC's for entire top electrode just:
    V_topBC(1:Nx+2, 1:Ny+2) = Va/Vt;
    
    
    %     V_topBC(ceil(num_cell/2):ceil(num_cell/2)+1, ceil(num_cell/2):ceil(num_cell/2)+1) = Va/Vt;  %only pts corresponding to tip are non-zero (rest BC are 0 EVEN THOUGH
    %ACTUAL V VALUE  MIGHT NOT BE 0, THERE IS NO BC ADDED TO RHS FOR NEUMANN BC'S.
    %other values already set to 0 earlier.
    
    iter = 1;
    error_np =  1;
    % Solver loop
    while(error_np > tolerance)
        
        %% Poisson Solve
        bV = SetbV_3D(p, epsilon);
        
        %solve for V
        oldV = V;
        
        %             newV = U\(L\bV);  %much faster to solve pre-factorized matrix. Not applicable to cont. eqn. b/c matrices keep changing.
        %          newV = AV\bV;
        [newV,~] = bicgstab(AV,bV, 10^-14, 1000, [], [], V); %SPECIFYING AN INITIAL GUESS IS CRUCIAL FOR SPEED!
        
        %NOTE: with scaling of AV and bV by 18 (diagonal amount), this
        %bicgstab converges much better!
        
        %NOTE: USING bicgstab as default (w/o specifying a tolerance), results in BAD results!!
        
        if(iter >1)
            V = newV*w + oldV*(1.-w);
        else
            V = newV;  %no mixing for 1st iter
        end
        
        %reshape the V vector into the V matrix
        
        V_matrix = reshape(V,Nz+1,Ny+1,Nx+1); %it is this ordering b/c V has order of z,x,y
        %need to permute the matrix to go back to x,y,z format
        V_matrix = permute(V_matrix, [3 2 1]);
        %---------------------------------------------------------------------------------
        %add on the BC's to get full potential matrix
        fullV(2:Nx+2,2:Ny+2, 2:Nz+2) = V_matrix;
        fullV(1:Nx+2,1:Ny+2,1) = V_bottomBC(1:Nx+2,1:Ny+2);
        fullV(1,2:Ny+2,2:Nz+2) = V_matrix(Nx+1,:,:);  %x bc's
        fullV(2:Nx+2,1,2:Nz+2) = V_matrix(:,Ny+1,:);
        %fill edges
        fullV(1,1,2:Nz+2) = V_matrix(Nx+1,1,:);  %left = right
         fullV(1,Ny+2,2:Nz+2) = V_matrix(1,Ny+1,:);  %this one seems
%         unneeded
        
        
        
%         %% Update net generation rate
% %         if(Va_cnt > 0)
% %             Up = 0;   %these only  include insides since want ctoo be consistent with i,j matrix indices
% %             %Up= zeros(num_elements,num_elements);
% % %             Un= Up;
% %         end

        
        %% Continuity equations solve
        Bernoulli_p_values = Calculate_Bernoullis_p(fullV);  %the values are returned as a struct
%         Bernoulli_n_values = Calculate_Bernoullis_n(fullV);
        Ap = SetAp_3D(p_mob_avged, Bernoulli_p_values);           
%         An = SetAn_3D(n_mob_avged, Bernoulli_n_values);
        bp = Setbp_3D(Bernoulli_p_values, p_mob, Up);
%         bn = Setbn_3D(Bernoulli_n_values, n_mob, Un);
        
        
        %NOTE TO SEE THE WHOLE SPARSE MATRICES: use :
        % full([name of sparse matrix])
        
        oldp = p;
        
        
%                   newp = Ap\bp;

        [newp, ~] = bicgstab(Ap, bp, 10^-14, 1000, [], [], p);  %Note: if don't specify an initial guess, bicgstab FAILS! 
        
   

%-------------------------------
%NOTE: BE VERY CAREFUL WITH  BICGSTAB, CHECK THE OUTPUT INFOR--> IT
%SOMETIMES FAILS TO CONVERGE AND GIVES BAD RESULTS WITHOUT GIVING A
%WARNING!
%---------------------------------
        
%         oldn = n;
        %          newn = An\bn;
%         [newn, ~] = bicgstab(An,bn, 10^-14, 1000);  %NOTE: using a lower tolerance, i..e 10^-14, makes this faster--> since more accuracy makes overall convergence of DD model faster
        %Note: seems 10^-14, is about the best level it can converge to
        
        % if get negative p's or n's, make them equal 0
        for i = 1:num_elements
            if(newp(i) <0.0)
                newp(i) = 0;
            end
%             if(newn(i) <0.0)
%                 newn(i) = 0;
%             end
        end
        
        old_error =  error_np;
        count = 0;
        error_np_matrix = zeros(1,num_elements); %need to reset the matrix b/c when more newp's become 0, need to remove that error element from matrix.
        for i = 1:num_elements  %recall that newp, oldp are stored as vectors
            if(newp(i) ~=0)
                count = count+1;  %counts number of non zero error calculations
                error_np_matrix(count) = (abs(newp(i)-oldp(i)))./abs(oldp(i));  %need the dot slash, otherwise it tries to do matrix operation! %ERROR SHOULD BE CALCULATED BEFORE WEIGHTING
            end
        end
        error_np = max(error_np_matrix)
        
        %auto decrease w if not converging
        if(error_np>= old_error)
            not_cnv_cnt = not_cnv_cnt+1;
        end
        if(not_cnv_cnt>100)
            w = w/2.;
            tolerance = tolerance*10;
            not_cnv_cnt = 0;  %reset the count
        end
        
        %w
        %tolerance
        
        %weighting
        p = newp*w + oldp*(1.-w);
%         n = newn*w + oldn*(1.-w);
        
        %% Apply continuity eqn BCs
        
        %reshape the vectors into matrices
        p_matrix = reshape(p,Nz+1,Ny+1,Nx+1);
        p_matrix = permute(p_matrix, [3 2 1]);
%         n_matrix = reshape(n,N,N,N);
      
        %-------------------------------------------------------------------------------------------------
        
        %add on the BC's to get full potential matrix
        fullp(2:Nx+2,2:Ny+2, 2:Nz+2) = p_matrix;
        fullp(1:Nx+2,1:Ny+2,1) = p_bottomBC(1:Nx+2,1:Ny+2);

        fullp(1,2:Ny+2,2:Nz+2) = p_matrix(Nx+1,:,:);  %x bc's
        fullp(2:Nx+2,1,2:Nz+2) = p_matrix(:,Ny+1,:);
        %fill edges
        fullp(1,1,2:Nz+2) = p_matrix(Nx+1,1,:);  %left = right
        fullp(1,Ny+2,2:Nz+2) = p_matrix(1,Ny+1,:);
        
        %add on the BC's to get full potential matrix
%         fulln(2:N+1,2:N+1, 2:N+1) = n_matrix;
%         fulln(1:N+2,1:N+2,1) = n_bottomBC;
%         fulln(1:N+2,1:N+2,N+2) = n_topBC;
%         fulln(1,2:N+1,2:N+1) = n_leftBC_x;
%         fulln(N+2,2:N+1,2:N+1) = n_rightBC_x;
%         fulln(2:N+1,1,2:N+1) = n_leftBC_y;
%         fulln(2:N+1,N+2,2:N+1) = n_rightBC_y;
%         %fill edges
%         fulln(1,1,2:N+1) = n_leftBC_x(1,:);
%         fulln(1,N+2,2:N+1) = n_leftBC_x(N,:);
%         fulln(N+2,1,2:N+1) = n_rightBC_x(1,:);
%         fulln(N+2,N+2,2:N+1) = n_rightBC_x(N,:);

        iter = iter +1;   
        
    
     
    end 
    
 
    % Calculate drift diffusion currents
    % Use the SG definition
    Bp_posX = Bernoulli_p_values.Bp_posX;
    Bp_negX = Bernoulli_p_values.Bp_negX;
    Bp_posY = Bernoulli_p_values.Bp_posY;
    Bp_negY = Bernoulli_p_values.Bp_negY;
    Bp_posZ = Bernoulli_p_values.Bp_posZ;
    Bp_negZ = Bernoulli_p_values.Bp_negZ;
    
%     Bn_posX = Bernoulli_n_values.Bn_posX;
%     Bn_negX = Bernoulli_n_values.Bn_negX;
%     Bn_posY = Bernoulli_n_values.Bn_posY;
%     Bn_negY = Bernoulli_n_values.Bn_negY;
%     Bn_posZ = Bernoulli_n_values.Bn_posZ;
%     Bn_negZ = Bernoulli_n_values.Bn_negZ;
    
    %the J(i+1,j+1) is to define J's (which are defined at mid cells, as rounded up integer)
 for i = 1:num_cell_x-1
        for j = 1:num_cell_y-1
            for k = 1:num_cell_z-1
                Jp_Z(i+1,j+1,k+1) = -(q*Vt*N_dos*mobil/dz)*(p_mob(i+1,j+1,k+1)*fullp(i+1,j+1,k+1)*Bp_negZ(i+1,j+1,k+1)-p_mob(i+1,j+1,k+1)*fullp(i+1,j+1,k)*Bp_posZ(i+1,j+1,k+1));
%                 Jn_Z(i+1,j+1,k+1) =  (q*Vt*N_dos*mobil/dz)*(n_mob(i+1,j+1,k+1)*fulln(i+1,j+1,k+1)*Bn_posZ(i+1,j+1,k+1)-n_mob(i+1,j+1,k+1)*fulln(i+1,j+1,k)*Bn_negZ(i+1,j+1,k+1));
                
                Jp_X(i+1,j+1,k+1) = -(q*Vt*N_dos*mobil/dx)*(p_mob(i+1,j+1,k+1)*fullp(i+1,j+1,k+1)*Bp_negX(i+1,j+1,k+1)-p_mob(i+1,j+1,k+1)*fullp(i,j+1,k+1)*Bp_posX(i+1,j+1,k+1));
%                 Jn_X(i+1,j+1,k+1) =  (q*Vt*N_dos*mobil/dx)*(n_mob(i+1,j+1,k+1)*fulln(i+1,j+1,k+1)*Bn_posX(i+1,j+1,k+1)-n_mob(i+1,j+1,k+1)*fulln(i,j+1,k+1)*Bn_negX(i+1,j+1,k+1));
                
                Jp_Y(i+1,j+1,k+1) = -(q*Vt*N_dos*mobil/dy)*(p_mob(i+1,j+1,k+1)*fullp(i+1,j+1,k+1)*Bp_negY(i+1,j+1,k+1)-p_mob(i+1,j+1,k+1)*fullp(i+1,j,k+1)*Bp_posY(i+1,j+1,k+1));
%                 Jn_Y(i+1,j+1,k+1) =  (q*Vt*N_dos*mobil/dy)*(n_mob(i+1,j+1,k+1)*fulln(i+1,j+1,k+1)*Bn_posY(i+1,j+1,k+1)-n_mob(i+1,j+1,k+1)*fulln(i+1,j,k+1)*Bn_negY(i+1,j+1,k+1));
            end        
        end
    end
    J_total_Z = Jp_Z;
    J_total_X = Jp_X;
    J_total_Y = Jp_Y;
    
    
    %Setup for JV curve
  if(Va_cnt>0)
        V_values(Va_cnt,:) = Va;
        J_total_Z_middle(Va_cnt) = J_total_Z(floor(num_cell_x/2),floor(num_cell_y/2),floor(num_cell_z/2));  %just store the J (in perpendicular to electrodes direction) in middle for the JV curve output
    end
    
    %Save data
    str = sprintf('%.2f',Va);
    filename = [str 'V.txt']
    fid = fopen(fullfile(filename),'w'); %fullfile allows to make filename from parts
    
    Up_matrix = reshape(Up,Nz+1,Ny+1,Nx+1);  %coming from z, y, x format
    Up_matrix = permute(Up_matrix, [3 2 1]);
%     Un_matrix = reshape(Un,N,N);
    
    if(Va_cnt > 0)
        for k = 2:num_cell_z
            i = floor(num_cell_x/2);     %JUST OUTPUT LINE PROFILE ALONG Z --> otherwise hard to compare results to 1D model
            j = floor(num_cell_y/2);
            fprintf(fid,'%.2e %.2e %.2e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\r\n', (i-1)*dx, (j-1)*dy, (k-1)*dz, Vt*fullV(i,j,k), N_dos*fullp(i,j,k), J_total_X(i,j,k), J_total_Z(i,j,k), Up_matrix(i-1,j-1,k-1), w, tolerance);
        end
    end
    fclose(fid);
    
    %% JV setup:
    file2 = fopen(fullfile('JV.txt'),'w');
    for i = 1:Va_cnt
        fprintf(file2, '%.8e %.8e \r\n', V_values(i,1), J_total_Z_middle(i));
    end
    fclose(file2);
    
    %% Final Plots: done for the last Va
    str = sprintf('%.2g', Va);
    
    Va_loop_time = toc(startLoop)  %time each Va
    
end

Program_time = toc(startProgram)

xslice = [  ceil((num_cell_x+1)/2)+1]  %note: am ensuring that including the mesh pots where tip occurs
yslice = [   ceil((num_cell_y+1)/2)+1]
zslice = [floor(num_cell_z/2), num_cell_z-1, num_cell_z, num_cell_z+1];
slice(fullV, xslice,yslice,zslice);

%to see just z planes
figure
xslice = [  ]  %note: am ensuring that including the mesh pots where tip occurs
yslice = [   ]
zslice = [2 3 floor(num_cell_z/2), num_cell_z-1, num_cell_z, num_cell_z+1];
slice(fullV, xslice,yslice,zslice);

%plot hole density
% figure
% xslice = [ floor(0.75*num_cell), floor(0.9*num_cell)]  %note: am ensuring that including the mesh pots where tip occurs
% yslice = [  floor(0.75*num_cell), floor(0.9*num_cell)]
% zslice = [floor(num_cell/2), num_cell-1, num_cell, num_cell+1];
% slice(fullp, xslice,yslice,zslice);

figure
xslice = [ ceil((num_cell_x+1)/2)+1]  %this puts the slices right at center of tip region
yslice = [  ceil((num_cell_y+1)/2)+1]
zslice = [floor(num_cell_z/2), num_cell_z-1, num_cell_z, num_cell_z];
slice(N_dos*fullp(:,:,1:num_cell_z), xslice,yslice,zslice);

%z holes slices only
figure
xslice = []  %this puts the slices right at center of tip region
yslice = [  ]
zslice = [2 3 floor(num_cell_z/2), num_cell_z-1, num_cell_z];
slice(N_dos*fullp(:,:,1:num_cell_z), xslice,yslice,zslice);  %plot NEGATIVE of J_total_Z just to get positive currents...

%z current slices only
figure
xslice = []  %this puts the slices right at center of tip region
yslice = [  ]
zslice = [2 3 floor(num_cell_z/2), num_cell_z-1, num_cell_z];
slice(-J_total_Z, xslice,yslice,zslice);  %plot NEGATIVE of J_total_Z just to get positive currents...


figure
xslice = [ ceil((num_cell_x+1)/2)+1]  %this puts the slices right at center of tip region
yslice = [  ceil((num_cell_y+1)/2)+1]
zslice = [floor(num_cell_z/2), num_cell_z-1, num_cell_z, num_cell_z];
slice(-J_total_Z, xslice,yslice,zslice);  %plot NEGATIVE of J_total_Z just to get positive currents...

%test the total Z current in each plane
I_Z = zeros(num_cell_z,1);
for k = 2:num_cell_z-1
    for i = 2:num_cell_x-1
        for j = 2:num_cell_y-1
             I_Z(k) = I_Z(k) + J_total_Z(i,j,k);
        end
    end
end


%plot p
% surf(1:N+2,1:N+2, log(N_dos*fullp))
% hold on
% surf(1:N+2,1:N+2, log(N_dos*fulln))
% hold off
% xlabel('Position ($m$)','interpreter','latex','FontSize',14);
% zlabel({'Log of carrier densities ($1/m^3$)'},'interpreter','latex','FontSize',14);

%plot line profiles of charge densities in the thickness direction
figure
plot(log(N_dos*fullp(1,1:Nz+2)))
hold on
% plot(log(N_dos*fulln(1,1:N+2)))
hold off
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Log of carrier densities ($1/m^3$)'},'interpreter','latex','FontSize',14);

figure
plot(Vt*fullV(1,1:Nz+2))
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Electric Potential (V)'},'interpreter','latex','FontSize',14);

