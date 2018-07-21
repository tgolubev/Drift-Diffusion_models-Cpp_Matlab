% Setup of 3D Poisson matrix using loop and spdiag

function AV = SetAV_3D(epsilon_avged)

global num_elements Nx Ny Nz dx dy dz

AV_val = zeros(2*num_elements, 11);   %this is a matrix which will just store the non-zero diagonals of 3D Poisson matrix

eps_X_avg = (dz^2/dx^2)*epsilon_avged.eps_X_avg./18.;  %scale
eps_Y_avg = (dz^2/dy^2)*epsilon_avged.eps_Y_avg./18.;  %to take into account possible dz dx dy difference, multipy by coefficient...
eps_Z_avg = epsilon_avged.eps_Z_avg./18.;

%--------------------------------------------------------------------------
%Lowest diagonal: X's Left PBC
values_cut = zeros(Nx+1, Ny+1, Nz+1);     %make i and j have N+1 b/c we are including the right bndry (n+1th) pt... for PBCs. Note: k still has N pts, b/c no pbc there
values_cut(1,1:Ny+1,1:Nz+1) = -eps_X_avg(2,2:Ny+1+1,2:Nz+1+1); %just index = 1 for x values--> since all correspond to left X pbc

%for dirhlet BC --> all elements corresponding to Nz+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);  
AV_val(1:num_elements, 1) = values_cut(:);  %recall lower diag's fill from the top, so don't need to worry about extra elements etc...


%--------------------------------------------------------------------------

%Lower diagonal:  X's
values_cut = zeros(Nx+1, Ny+1, Nz+1);     %make i and j have N+1 b/c we are including the right bndry (n+1th) pt... for PBCs. Note: k still has N pts, b/c no pbc there
values_cut(1:Nx+1-1,1:Ny+1,1:Nz+1) = -eps_X_avg(2+1:Nx+1+1,2:Ny+1+1,2:Nz+1+1); 

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);  %this changes order of indices, so have z, y, x --> z are the rows--> then when use (:), will have correct order in AV_val
AV_val(1:num_elements, 2) = values_cut(:); 


%--------------------------------------------------------------------------
%Lower diagonal: Y's LEFT PBC's
values_cut = zeros(Nx+1, Ny+1, Nz+1);      %RESET THE SIZE OF VALUES_CUT: this is
                                      %important!!                                      
  %note: made the j of N+1 size-->  so elements unfilled --> this corresponds to the 0's subblocks
values_cut(1:Nx+1,1,1:Nz+1) = -eps_Y_avg(2:Nx+1+1, 2, 2:Nz+1+1);  %Note: made y values be just 1 index --> b/c this pbc--> is just always refering to j = 1 left pbc.

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
AV_val(1:num_elements, 3) = values_cut(:);  %is N^3 lenght, but the last block will just have 0's....

%--------------------------------------------------------------------------
%Lower diagonal: Y's
values_cut = zeros(Nx+1, Ny+1, Nz+1);          %RESET THE SIZE OF VALUES_CUT: this is
                                      %important!!                                    
  %note: made the j of N+1 size-->  so when leave the last Nth elements,
  %unfilled--> this corresponds to the 0's subblocks

values_cut(1:Nx+1,1:Ny+1-1,1:Nz+1) = -eps_Y_avg(2:Nx+1+1,2+1:Ny+1+1,2:Nz+1+1);  


%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
AV_val(1:num_elements, 4) = values_cut(:);  %is N^3 lenght, but the last block will just have 0's....

%--------------------------------------------------------------------------
%main lower diagonal (below main diagonal)
values_cut = zeros(Nx+1, Ny+1, Nz+1); 
values_cut(1:Nx+1,1:Ny+1,1:Nz-1+1) = -eps_Z_avg(2:Nx+1+1,2:Ny+1+1,2+1:Nz+1+1); 

%FOR DIRICLET BC'S, for all i,j, and k = Nth index, need 0's
values_cut(:,:,Nz) = 0;

values_cut = permute(values_cut, [3 2 1]);
AV_val(1:num_elements, 5) = values_cut(:); 


%--------------------------------------------------------------------------
%main diagonal
values_cut = zeros(Nx+1, Ny+1, Nz+1); 
values_cut(1:Nx+1,1:Ny+1,1:Nz+1) = eps_X_avg(2:Nx+1+1, 2:Ny+1+1,2:Nz+1+1) + eps_X_avg(2+1:Nx+1+1+1, 2:Ny+1+1,2:Nz+1+1) + eps_Y_avg(2:Nx+1+1, 2:Ny+1+1,2:Nz+1+1) + eps_Y_avg(2:Nx+1+1,2+1:Ny+1+1+1,2:Nz+1+1) + ...
                              +eps_Z_avg(2:Nx+1+1, 2:Ny+1+1,2:Nz+1+1) + eps_Z_avg(2:Nx+1+1,2:Ny+1+1,2+1:Nz+1+1+1);  %note: k values go only to N+1 or max to N+1+1...., 1 less than X and Y

%FOR DIRICLET BC'S, for all i,j, and k = N+1th index, need 1's
values_cut(:,:,Nz+1) = 1;
                          
values_cut = permute(values_cut, [3 2 1]);
AV_val(1:num_elements,6) = values_cut(:);

%--------------------------------------------------------------------------
%main uppper diagonal
values_cut = zeros(Nx+1, Ny+1, Nz+1);
values_cut(1:Nx+1, 1:Ny+1,1:Nz-1+1) = -eps_Z_avg(2:Nx+1+1, 2:Ny+1+1, 2+1:Nz+1+1);  %2+1 b/c is k+1+1

%don't need anything here for Dirichlet BC's since those elements are
%already 0

values_cut = permute(values_cut, [3 2 1]);
AV_val(2:num_elements+1,7) = values_cut(:);  %+1 here b/c of way spdiags fills


%--------------------------------------------------------------------------
%upper diagonal   Y's
values_cut = zeros(Nx+1, Ny+1, Nz+1);
values_cut(1:Nx+1, 1:Ny+1-1,1:Nz+1) = -eps_Y_avg(2:Nx+1+1, 2+1:Ny+1+1, 2:Nz+1+1);  %note is j+1+1

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
AV_val(1+(Nz+1):num_elements+(Nz+1), 8) = values_cut(:);  %must shift by N+1


%--------------------------------------------------------------------------
%upper diagonal: right Y PBC's
values_cut = zeros(Nx+1, Ny+1, Nz+1);
values_cut(1:Nx+1, 1, 1:Nz+1) = -eps_Y_avg(2:Nx+1+1, Ny+1+1, 2:Nz+1+1);  %j = N+1+1 b/c corresponding to right pbc's

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
AV_val(1+(Nz+1)*Ny:num_elements+(Nz+1)*Ny, 9) = values_cut(:);  %must shift by N+1

%--------------------------------------------------------------------------
%upper diagonal: X's
values_cut = zeros(Nx+1, Ny+1, Nz+1);
values_cut(1:Nx+1-1, 1:Ny+1,1:Nz+1) = -eps_X_avg(2+1:Nx+1+1, 2:Ny+1+1, 2:Nz+1+1); %i+1+1

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
 values_cut(:,:,Nz+1) = 0; %this is correct
 
values_cut = permute(values_cut, [3 2 1]);
AV_val(1+(Nz+1)*(Ny+1):num_elements + (Nz+1)*(Ny+1) ,10) = values_cut(:);  %shifted by N*N+1 b/c N z values and N+1 y values to shift over by.

%--------------------------------------------------------------------------
%far upper diagonal: X right PBC's
values_cut = zeros(Nx+1, Ny+1, Nz+1);
values_cut(1, 1:Ny+1,1:Nz+1) = -eps_X_avg(Nx+1+1, 2:Ny+1+1, 2:Nz+1+1); %i = N+1+1 b/c always corresonding to far right BC

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
AV_val(1+(Nx)*(Nz+1)*(Ny+1):num_elements + (Nx)*(Nz+1)*(Ny+1) ,11) = values_cut(:);  %shifted b/c spdiags fills upper diags from the bottom

%--------------------------------------------------------------------------
AV = spdiags(AV_val, [-(Nx)*(Nz+1)*(Ny+1) -(Nz+1)*(Ny+1) -Ny*(Nz+1) -(Nz+1) -1 0 1 Nz+1 Ny*(Nz+1) (Nz+1)*(Ny+1) (Nx)*(Nz+1)*(Ny+1)], num_elements, num_elements); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.


