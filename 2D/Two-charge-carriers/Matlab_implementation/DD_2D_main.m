%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         2D Drift Diffusion for Electrons and Holes with Finite Differences
%
%             -Scharfetter-Gummel discretization
%             -decoupled Gummel iterations method
%
%              Created by: Timofey Golubev (2018.05.29)
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
%
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
global num_cell N num_elements Vt N_dos p_topBC p_leftBC p_rightBC p_bottomBC Cn CV Cp n_leftBC  n_rightBC
global n_bottomBC n_topBC V_leftBC V_bottomBC V_topBC V_rightBC G_max

%% Physical Constants

q =  1.60217646*10^-19;         %elementary charge (C)
kb = 1.3806503*10^-23;          %Boltzmann const. (J/k)
T = 296.;                       %temperature (K)
epsilon_0 =  8.85418782*10^-12; %F/m
Vt = (kb*T)/q;

%% Simulation Setup

%Voltage sweep loop
Va_min = -0.5;            %volts
Va_max = 1.0;
increment = 0.01;         %by which to increase V
num_V = floor((Va_max-Va_min)/increment)+1;   %number of V points

%Simulation parameters
w_eq = 0.05;               %linear mixing factor for 1st convergence (0 applied voltage, no generation equilibrium case)
w_i = 0.2;                 %starting linear mixing factor for Va_min (will be auto decreased if convergence not reached)
tolerance = 5*10^-12;        %error tolerance
tolerance_i =  5*10^-12;     %initial error tolerance, will be increased if can't converge to this level

%% System Setup

L = 10.0000001e-9;     %there's some integer rounding issue, so use this .0000001
dx = 1e-9;                        %mesh size
num_cell = floor(L/dx);
N = num_cell -1;       %number of INTERIOR mesh points (total mesh pts = num_cell +1 b/c matlab indixes from 1)
num_elements = N^2;  %NOTE: this will specify number of elements in the solution vector V which = (num_cell +1 -2)^2 b/c we are in 2D

%Electronic density of states of holes and electrons
N_VB = 10^24;         %density of states in valence band (holes)
N_CB = 10^24;         %density of states in conduction bands (electrons)
E_gap = 1.5;          %bandgap of the active layer(in eV)
N_dos = 10^24.;       %scaling factor helps CV be on order of 1

%injection barriers
inj_a = 0.2;	%at anode
inj_c = 0.1;	%at cathode

%work functions of anode and cathode
WF_anode = 4.8;
WF_cathode = 3.7;

Vbi = WF_anode - WF_cathode +inj_a +inj_c;  %built-in field

G_max = 4*10^27;

%% Define matrices of system parameters

% Relative dielectric constant matrix (can be position dependent)
%Epsilons are defined at 1/2 integer points, so epsilons inside
%the cells, not at cell boundaries
%will use indexing: i + 1/2 is defined as i+1 for the index
epsilon = zeros(num_cell+1, num_cell +1);
for i = 1:num_cell+1
    epsilon(i,:) = 3.0;
end

% Define mobilities matrix (can be position dependent)
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

%Scaling coefficients
Cp = dx^2/(Vt*N_dos*mobil);          %note: scaled p_mob and n_mob are inside matrices
Cn = dx^2/(Vt*N_dos*mobil);
CV = N_dos*dx^2*q/(epsilon_0*Vt);    %relative permitivity was moved into the matrix

%% Define Poisson equation boundary conditions and initial conditions

% Initial conditions
V_bottomBC = -((Vbi)/(2*Vt)-inj_a/Vt);
V_topBC = (Vbi)/(2*Vt)-inj_c/Vt;
diff = (V_topBC - V_bottomBC)/num_cell;
V(1:N) = V_bottomBC + diff;  %define V's corresponding to 1st subblock here (1st interior row of system)
index = N;
for j = 2:N  %corresponds to z coord
    index = index +1;
    V(index) = V(index-1) + diff;
    for i = 2:N  %elements along the x direction assumed to have same V
        index = index +1;
        V(index) = V(index-1);
    end
end

%reshape the V vector into the V matrix
V_matrix = reshape(V,N,N);

%Side BCs (insulating)
for j = 1:N
    V_leftBC(j) = V((j-1)*N + 1);  %just corresponds to values of V in the 1st subblock
    V_rightBC(j) = V(j*N);
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

%form matrices for easy filling of bp
n_matrix = reshape(n,N,N);
p_matrix = reshape(p,N,N);

%side boundary conditions %insulating bc's,: charge at bndry = charge just inside
for j = 1:N  
    p_leftBC(j) = p((j-1)*N + 1);
    p_rightBC(j)= p(j*N);
    n_leftBC(j) = n((j-1)*N + 1);
    n_rightBC(j)= n(j*N);
end

% Set up Poisson matrix equation
AV = SetAV_2D(epsilon);
%spy(AV);  %allows to see matrix structure, very useful!

%% Main voltage loop
Va_cnt = 0;
for Va_cnt = 0:num_V +1
    not_converged = false;
    not_cnv_cnt = 0;
    
    %stop the calculation if tolerance becomes too high
    if(tolerance >10^-5)
        break
    end
    
    %1st iteration is to find the equilibrium values (no generation rate)
    if(Va_cnt ==0)
        tolerance = tolerance*10^2;       %relax tolerance for equil convergence
        w = w_eq;                         %use smaller mixing factor for equil convergence
        Va = 0;
        Up = zeros(num_cell+1,num_cell+1);
        Un= Up;
    end
    if(Va_cnt ==1)
        tolerance = tolerance_i;       %reset tolerance back
        w=w_i;
    end
    if(Va_cnt >0)
        Va = Va_min+increment*(Va_cnt-1);     %set Va value
        %Va = Va_max-increment*(Va_cnt-1);    %decrease Va by increment in each iteration
    end
    
    %Voltage boundary conditions
    V_bottomBC = -((Vbi  -Va)/(2*Vt)-inj_a/Vt);
    V_topBC = (Vbi- Va)/(2*Vt) - inj_c/Vt;
      
    iter = 1;
    error_np =  1;
    % Solver loop
    while(error_np > tolerance)
        %% Poisson Solve
        bV = SetbV_2D(p_matrix, n_matrix, epsilon);
        
        %solve for V
        oldV = V;
        newV = AV\bV;
        
        if(iter >1)
            V = newV*w + oldV*(1.-w);
        else
            V = newV;  %no mixing for 1st iter
        end
        
        %reshape the V vector into the V matrix
        V_matrix = reshape(V,N,N);
        
        %Update side boundary conditions
        for j = 1:N
            V_leftBC(j) = V((j-1)*N + 1); 
            V_rightBC(j) = V(j*N);
        end
        
        %add on the BC's to get full potential matrix
        fullV(2:N+1,2:N+1) = V_matrix;
        fullV(1:N+2,1) = V_bottomBC;
        fullV(:,N+2) = V_topBC;
        fullV(1,2:N+1) = V_leftBC;
        fullV(N+2,2:N+1) = V_rightBC;
        fullV(N+2, N+2) = V_topBC;    %far exterior corner has same V as rest of the top
        
        %% Update net generation rate
        if(Va_cnt > 0)
            G = GenerationRate();
            Up(1:N,1:N) = G(1:N,1:N);  %only defined on interior of the device
            %Up= zeros(num_elements,num_elements);
            Un= Up;
        end
        
        %% Continuity equations solve
        Bernoulli_p_values = Calculate_Bernoullis_p(fullV);  %the values are returned as a struct
        Bernoulli_n_values = Calculate_Bernoullis_n(fullV);
        Ap = SetAp_2D(p_mob, Bernoulli_p_values);           
        An = SetAn_2D(n_mob, Bernoulli_n_values);
        bp = Setbp_2D(Bernoulli_p_values, p_mob, Up);
        bn = Setbn_2D(Bernoulli_n_values, n_mob, Un);
        
        %NOTE TO SEE THE WHOLE SPARSE MATRICES: use :
        % full([name of sparse matrix])
        
        oldp = p;
        newp = Ap\bp;
        
        oldn = n;
        newn = An\bn;
        
        % if get negative p's or n's, make them equal 0
        for i = 1:num_elements
            if(newp(i) <0.0)
                newp(i) = 0;
            end
            if(newn(i) <0.0)
                newn(i) = 0;
            end
        end
        
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
        
        %auto decrease w if not converging
        if(error_np>= old_error)
            not_cnv_cnt = not_cnv_cnt+1;
        end
        if(not_cnv_cnt>2000)
            w = w/2.;
            tolerance = tolerance*10;
            not_cnv_cnt = 0;  %reset the count
        end
        
        w
        tolerance
        
        %weighting
        p = newp*w + oldp*(1.-w);
        n = newn*w + oldn*(1.-w);
        
        %% Apply continuity eqn BCs
        
        %reshape the vectors into matrices
        p_matrix = reshape(p,N,N);
        n_matrix = reshape(n,N,N);
        
        %Update side boundary conditions
        for j = 1:N
            p_leftBC(j) = p((j-1)*N + 1); %just corresponds to values of V in the 1st subblock
            p_rightBC(j)= p(j*N);
        end
        
        %add on the BC's to get full potential matrix
        fullp(2:N+1,2:N+1) = p_matrix;
        fullp(1:N+2,1) = p_bottomBC;
        fullp(:,N+2) = p_topBC;
        fullp(1,2:N+1) = p_leftBC;
        fullp(N+2,2:N+1) = p_rightBC;
        fullp(N+2, N+2) = p_topBC;  
        
        for j = 1:N
            n_leftBC(j) = n_matrix(N,j);
            n_rightBC(j)= n_matrix(1,j);
        end
        
        %add on the BC's to get full potential matrix
        fulln(2:N+1,2:N+1) = n_matrix;
        fulln(1:N+2,1) = n_bottomBC;
        fulln(:,N+2) = n_topBC;
        fulln(1,2:N+1) =n_leftBC;
        fulln(N+2,2:N+1) = n_rightBC;
        fulln(N+2, N+2) = n_topBC;  
        
        iter = iter +1   
    end
    
    % Calculate drift diffusion currents
    % Use the SG definition
    Bp_posX = Bernoulli_p_values.Bp_posX;
    Bp_negX = Bernoulli_p_values.Bp_negX;
    Bp_posZ = Bernoulli_p_values.Bp_posZ;
    Bp_negZ = Bernoulli_p_values.Bp_negZ;
    Bn_posX = Bernoulli_n_values.Bn_posX;
    Bn_negX = Bernoulli_n_values.Bn_negX;
    Bn_posZ = Bernoulli_n_values.Bn_posZ;
    Bn_negZ = Bernoulli_n_values.Bn_negZ;
    
    %the J(i+1,j+1) is to define J's (which are defined at mid cells, as rounded up integer)
    for i = 1:num_cell-1
        for j = 1:num_cell-1
            Jp_Z(i+1,j+1) = -(q*Vt*N_dos*mobil/dx)*(p_mob(i+1,j+1)*fullp(i+1,j+1)*Bp_negZ(i+1,j+1)-p_mob(i+1,j+1)*fullp(i+1,j)*Bp_posZ(i+1,j+1));
            Jn_Z(i+1,j+1) =  (q*Vt*N_dos*mobil/dx)*(n_mob(i+1,j+1)*fulln(i+1,j+1)*Bn_posZ(i+1,j+1)-n_mob(i+1,j+1)*fulln(i+1,j)*Bn_negZ(i+1,j+1));
            
            Jp_X(i+1,j+1) = -(q*Vt*N_dos*mobil/dx)*(p_mob(i+1,j+1)*fullp(i+1,j+1)*Bp_negX(i+1,j+1)-p_mob(i+1,j+1)*fullp(i,j+1)*Bp_posX(i+1,j+1));
            Jn_X(i+1,j+1) =  (q*Vt*N_dos*mobil/dx)*(n_mob(i+1,j+1)*fulln(i+1,j+1)*Bn_posX(i+1,j+1)-n_mob(i+1,j+1)*fulln(i,j+1)*Bn_negX(i+1,j+1));
        end
    end
    J_total_Z = Jp_Z + Jn_Z;
    J_total_X = Jp_X + Jn_X;
    
    %Setup for JV curve
    if(Va_cnt>0)
        V_values(Va_cnt,:) = Va;
        J_total_Z_middle(Va_cnt) = J_total_Z(floor(num_cell/2),floor(num_cell/2));  %just store the J (in perpendicular to electrodes direction) in middle for the JV curve output
    end
    
    %Save data
    str = sprintf('%.2f',Va);
    filename = [str 'V.txt']
    fid = fopen(fullfile(filename),'w'); %fullfile allows to make filename from parts
    
    if(Va_cnt ==0)
        equil = fopen(fullfile('Equil.txt'),'w');
        for j = 2:num_cell
            i = floor(num_cell/2);
            fprintf(equil,'%.2e %.2e %.8e %.8e %.8e \r\n ',(i-1)*dx, (j-1)*dx, Vt*fullV(i,j), N_dos*fullp(i,j), N_dos*fulln(i,j));
        end
        fclose(equil);
    end
    
    if(Va_cnt > 0)
        for j = 2:num_cell
            i = floor(num_cell/2);     %JUST OUTPUT LINE PROFILE ALONG Z --> otherwise hard to compare results to 1D model
            %for j = 2:num_cell    %use if want to output values along the entire grid
            fprintf(fid,'%.2e %.2e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\r\n', (i-1)*dx, (j-1)*dx, Vt*fullV(i,j), N_dos*fullp(i,j), N_dos*fulln(i,j), J_total_X(i,j), J_total_Z(i,j), Up(i-1,j-1), w, tolerance);
            %end
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
    
end

%plot V
surf(1:N+2,1:N+2, Vt*fullV)
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Electric Potential (V)'},'interpreter','latex','FontSize',14);

figure
%plot p
surf(1:N+2,1:N+2, log(N_dos*fullp))
hold on
surf(1:N+2,1:N+2, log(N_dos*fulln))
hold off
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
zlabel({'Log of carrier densities ($1/m^3$)'},'interpreter','latex','FontSize',14);

%plot line profiles of charge densities in the thickness direction
figure
plot(log(N_dos*fullp(1,1:N+2)))
hold on
plot(log(N_dos*fulln(1,1:N+2)))
hold off
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Log of carrier densities ($1/m^3$)'},'interpreter','latex','FontSize',14);

figure
plot(Vt*fullV(1,1:N+2))
xlabel('Position ($m$)','interpreter','latex','FontSize',14);
ylabel({'Electric Potential (V)'},'interpreter','latex','FontSize',14);
