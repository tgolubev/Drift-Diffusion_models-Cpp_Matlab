% Setup of 3D hole continuity matrix using loop and spdiag

function Ap = SetAp_3D(p_mob_avged, Bernoulli_p_values)

global num_elements Nx Ny Nz dx dy dz

%extract variables from struct, for brevity in eqns
Bp_posX = Bernoulli_p_values.Bp_posX;
Bp_negX = Bernoulli_p_values.Bp_negX;
Bp_posY = Bernoulli_p_values.Bp_posY;
Bp_negY = Bernoulli_p_values.Bp_negY;
Bp_posZ = Bernoulli_p_values.Bp_posZ;
Bp_negZ = Bernoulli_p_values.Bp_negZ;

p_mob_X_avg = (dz^2/dx^2)*p_mob_avged.p_mob_X_avg;  %to take into account possible dz dx dy difference, multipy by coefficient...
p_mob_Y_avg = (dz^2/dy^2)*p_mob_avged.p_mob_Y_avg;
p_mob_Z_avg = p_mob_avged.p_mob_Z_avg;

%NOTE: ALL THE VALUES must start with index 2, b/c the Bernoulli's start
%with 2!!

Ap_val = zeros(2*num_elements, 11);   %this is a matrix which will just store the non-zero diagonals of 2D hole continuity eqn

%--------------------------------------------------------------------------
%Lowest diagonal: X's Left PBC
values = -p_mob_X_avg.*Bp_posX;
values_cut = zeros(Nx+1, Ny+1, Nz+1);     %make i and j have N+1 b/c we are including the right bndry (n+1th) pt... for PBCs. Note: k still has N pts, b/c no pbc there
values_cut(1,1:Ny+1,1:Nz+1) = values(2,2:Ny+1+1,2:Nz+1+1); %just index = 2 for x values--> since all correspond to left X pbc

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);  
Ap_val(1:num_elements, 1) = values_cut(:);  %recall lower diag's fill from the top, so don't need to worry about extra elements etc...


%--------------------------------------------------------------------------

%Lower diagonal:  X's
values = -p_mob_X_avg.*Bp_posX;
values_cut = zeros(Nx+1, Ny+1, Nz+1);     %make i and j have N+1 b/c we are including the right bndry (n+1th) pt... for PBCs. Note: k still has N pts, b/c no pbc there
values_cut(1:Nx+1-1,1:Ny+1,1:Nz+1) = values(2+1:Nx+1+1,2:Ny+1+1,2:Nz+1+1); 

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);  %this changes order of indices, so have z, y, x --> z are the rows--> then when use (:), will have correct order in Ap_val
Ap_val(1:num_elements, 2) = values_cut(:); 


%--------------------------------------------------------------------------
%Lower diagonal: Y's LEFT PBC's
values = -p_mob_Y_avg.*Bp_posY;
values_cut = zeros(Nx+1, Ny+1, Nz+1);      %RESET THE SIZE OF VALUES_CUT: this is
                                      %important!!
                            
  %note: made the j of N+1 size-->  so elements unfilled --> this corresponds to the 0's subblocks

values_cut(1:Nx+1,1,1:Nz+1) = values(2:Nx+1+1, 2, 2:Nz+1+1);  %Note: made y values be just 1 index --> b/c this pbc--> is just always refering to j = 1 left pbc.

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
Ap_val(1:num_elements, 3) = values_cut(:);  %is N^3 lenght, but the last block will just have 0's....

%--------------------------------------------------------------------------
%Lower diagonal: Y's
values = -p_mob_Y_avg.*Bp_posY;
values_cut = zeros(Nx+1, Ny+1, Nz+1);          %RESET THE SIZE OF VALUES_CUT: this is
                                      %important!!
                                      
  %note: made the j of N+1 size-->  so when leave the last Nth elements,
  %unfilled--> this corresponds to the 0's subblocks

values_cut(1:Nx+1,1:Ny+1-1,1:Nz+1) = values(2:Nx+1+1,2+1:Ny+1+1,2:Nz+1+1);  

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
Ap_val(1:num_elements, 4) = values_cut(:);  %is N^3 lenght, but the last block will just have 0's....

%--------------------------------------------------------------------------
%main lower diagonal (below main diagonal)
values = -p_mob_Z_avg.*Bp_posZ;
values_cut = zeros(Nx+1, Ny+1, Nz+1); 
values_cut(1:Nx+1,1:Ny+1,1:Nz-1+1) = values(2:Nx+1+1,2:Ny+1+1,2+1:Nz+1+1);


%FOR DIRICLET BC'S, for all i,j, and k = Nth index, need 0's
values_cut(:,:,Nz) = 0;

values_cut = permute(values_cut, [3 2 1]);
Ap_val(1:num_elements, 5) = values_cut(:); 

%--------------------------------------------------------------------------
%main diagonal
values1 = p_mob_Z_avg.*Bp_negZ + p_mob_Y_avg.*Bp_negY + p_mob_X_avg.*Bp_negX;
values2 = p_mob_X_avg.*Bp_posX;  %these have +1+1 in some elements, so need to be seperate
values3 = p_mob_Y_avg.*Bp_posY;
values4 = p_mob_Z_avg.*Bp_posZ;

values_cut = zeros(Nx+1, Ny+1, Nz+1); 
values_cut(1:Nx+1,1:Ny+1,1:Nz+1) = values1(2:Nx+1+1, 2:Ny+1+1,2:Nz+1+1) + values2(2+1:Nx+1+1+1, 2:Ny+1+1,2:Nz+1+1) + values3(2:Nx+1+1,2+1:Ny+1+1+1,2:Nz+1+1) + values4(2:Nx+1+1,2:Ny+1+1,2+1:Nz+1+1+1);


%FOR DIRICLET BC'S, for all i,j, and k = N+1th index, need 1's
values_cut(:,:,Nz+1) = 1;

values_cut = permute(values_cut, [3 2 1]);
Ap_val(1:num_elements,6) = values_cut(:);

%--------------------------------------------------------------------------
%main uppper diagonal
values = -p_mob_Z_avg.*Bp_negZ;
values_cut = zeros(Nx+1, Ny+1, Nz+1);
values_cut(1:Nx+1, 1:Ny+1,1:Nz-1+1) = values(2:Nx+1+1, 2:Ny+1+1, 2+1:Nz+1+1);  %2+1 b/c is k+1+1

%DON'T need to do anything here for diriclet Bc's, b/c subblock corners on
%this side are already 0

values_cut = permute(values_cut, [3 2 1]);
Ap_val(2:num_elements+1,7) = values_cut(:);  %+1 here b/c of way spdiags fills


%--------------------------------------------------------------------------
%upper diagonal   Y's
values = -p_mob_Y_avg.*Bp_negY;
values_cut = zeros(Nx+1, Ny+1, Nz+1);
values_cut(1:Nx+1, 1:Ny+1-1,1:Nz+1) = values(2:Nx+1+1, 2+1:Ny+1+1, 2:Nz+1+1);  %note is j+1+1

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
Ap_val(1+(Nz+1):num_elements+(Nz+1), 8) = values_cut(:);  %must shift by Nz+1 b/c N+1 elements in each subblock (due to including Neuman BC's have the +1)


%--------------------------------------------------------------------------
%upper diagonal: right Y PBC's
values = -p_mob_Y_avg.*Bp_negY;
values_cut = zeros(Nx+1, Ny+1, Nz+1);
values_cut(1:Nx+1, 1, 1:Nz+1) = values(2:Nx+1+1, Ny+1+1, 2:Nz+1+1);  %j = N+1+1 b/c corresponding to right pbc's

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
Ap_val(1+(Nz+1)*Ny:num_elements+(Nz+1)*Ny, 9) = values_cut(:);  %(Nz+1)*Ny gives distance from 1st element to the location of PBCs

%--------------------------------------------------------------------------
%upper diagonal: X's
values = -p_mob_X_avg.*Bp_negX;
values_cut = zeros(Nx+1, Ny+1, Nz+1);
values_cut(1:Nx+1-1, 1:Ny+1,1:Nz+1) = values(2+1:Nx+1+1, 2:Ny+1+1, 2:Nz+1+1); %i+1+1

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
Ap_val(1+(Nz+1)*(Ny+1):num_elements + (Nz+1)*(Ny+1) ,10) = values_cut(:);  %shifted by 1 whole subblock group (N+1)^2

%--------------------------------------------------------------------------
%far upper diagonal: X right PBC's
values = -p_mob_X_avg.*Bp_negX;
values_cut = zeros(Nx+1, Ny+1, Nz+1);
values_cut(1, 1:Ny+1,1:Nz+1) = values(Nx+1+1, 2:Ny+1+1, 2:Nz+1+1); %i = N+1+1 b/c always corresonding to far right BC

%for dirhlet BC --> all elements corresponding to N+1 = 0 --> b/c that's
%just  BC eqn
values_cut(:,:,Nz+1) = 0;

values_cut = permute(values_cut, [3 2 1]);
Ap_val(1+(Nx)*(Nz+1)*(Ny+1):num_elements + (Nx)*(Nz+1)*(Ny+1) ,11) = values_cut(:);  %shifted b/c spdiags fills upper diags from the bottom. Shifted by 1 less subblock group than total num elements

%--------------------------------------------------------------------------
Ap = spdiags(Ap_val, [-(Nx)*(Nz+1)*(Ny+1) -(Nz+1)*(Ny+1) -Ny*(Nz+1) -(Nz+1) -1 0 1 Nz+1 Ny*(Nz+1) (Nz+1)*(Ny+1) (Nx)*(Nz+1)*(Ny+1)], num_elements, num_elements); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.

