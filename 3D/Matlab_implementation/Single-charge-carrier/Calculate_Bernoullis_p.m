function Bernoulli_p_values = Calculate_Bernoullis_p(fullV)
global num_cell_x num_cell_y num_cell_z 

%This will return a vector/array containing valuess for all 4 Bernoulli p fncs as follows:
%Bernoulli_values = [Bp_posX  Bp_negX Bp_posZ Bp_negZ]

%matlab suggests to not preallocate these--> b/c are output of the
%function. Already allocated when call the function.
% Bp_posX = zeros(num_cell+1,num_cell+1);
% Bp_negX = zeros(num_cell+1,num_cell+1);
% Bp_posZ = zeros(num_cell+1,num_cell+1);
% Bp_negZ = zeros(num_cell+1,num_cell+1);
dV_X = zeros(num_cell_x+2,num_cell_y+2,num_cell_z+2);
dV_Y = zeros(num_cell_x+2,num_cell_y+2,num_cell_z+2);
dV_Z = zeros(num_cell_x+2,num_cell_y+2,num_cell_z+2);  %now need to add bc case for z's --> b/c including that bndry element
for k = 2:num_cell_z+1
    for j = 2:num_cell_y+1
        for i = 2:num_cell_x+1
            dV_X(i,j,k) = fullV(i,j,k) - fullV(i-1,j,k);   %NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
            dV_Y(i,j,k) = fullV(i,j,k) - fullV(i,j-1,k);
            dV_Z(i,j,k) = fullV(i,j,k) - fullV(i,j,k-1);
            
            %add boundary case for num_cell+2 on Z value--> that dV = 0
            %from the Neuman boundary condition! And where we do have tip,
            %we assume tip voltage is same at +1 into tip, as at tip
            %contact bndry pt.
            dV_X(i,j, num_cell_z+2) = 0;
            dV_Y(i,j, num_cell_z+2) = 0;
            dV_Z(i,j, num_cell_z+2) = 0;
            
            %add on the wrap around dV's
            dV_X(i,num_cell_y+2,k) = fullV(i,1,k) - fullV(i-1,1,k);   %NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
            dV_Y(i,num_cell_y+2,k) = fullV(i,1,k) - fullV(i,num_cell_y+1,k);
            dV_Z(i,num_cell_y+2,k) = fullV(i,1,k) - fullV(i,1,k-1);
            
        end
        %add on the wrap around dV's and right boundary pt dV's
        dV_X(num_cell_x+2,j,k) = fullV(1,j,k) - fullV(num_cell_x+1,j,k);   %NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
        dV_Y(num_cell_x+2,j,k) = fullV(1,j,k) - fullV(1,j-1,k);  %use 1's b/c num_cell+2 is same as the 1st pt--> PBC's
        dV_Z(num_cell_x+2,j,k) = fullV(1,j,k) - fullV(1,j,k-1);
        
        
    end
    dV_X(num_cell_x+2,num_cell_y+2,k) = fullV(1,j,k) - fullV(num_cell_x+1,j,k);   %NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
    dV_Y(num_cell_x+2,num_cell_y+2,k) = fullV(1,j,k) - fullV(1,j-1,k);  %use 1's b/c num_cell+2 is same as the 1st pt--> PBC's
    dV_Z(num_cell_x+2,num_cell_y+2,k) = fullV(1,j,k) - fullV(1,j,k-1);
    
    
end

% matrix operations for Bernoulli's, SAVES SIGNIFICANT CPU then overwrite the elements 
%which need overwriting...
%BE CAREFUL!, NEED TO USE THE . to indicate element-by-element operations
Bp_posX = dV_X./(exp(dV_X) - 1.0);
Bp_negX = Bp_posX.*exp(dV_X);

Bp_posY = dV_Y./(exp(dV_Y) - 1.0);
Bp_negY = Bp_posY.*exp(dV_Y);

Bp_posZ = dV_Z./(exp(dV_Z) - 1.0);
Bp_negZ = Bp_posZ.*exp(dV_Z);

% Bp_posX = ones(num_cell+1,num_cell+1,num_cell+1);
% Bp_negX = ones(num_cell+1,num_cell+1,num_cell+1);
% Bp_posY = ones(num_cell+1,num_cell+1,num_cell+1);
% Bp_negY = ones(num_cell+1,num_cell+1,num_cell+1);
% Bp_posZ = ones(num_cell+1,num_cell+1,num_cell+1);
% Bp_negZ = ones(num_cell+1,num_cell+1,num_cell+1);

%I THINK THIS IS STILL OK TO USE--> only a ffefcts the small dV's--> i.e.
%which are negligible. For larger dV's, we use the real Bernoulli eqn's...
for i = 2:num_cell_x+2  %num_cell +2 for i and j b/c of the PBC's/including boundary pt--> but are setting to 1 anyway..., so doesn't matter much
    for j = 2:num_cell_y+2
        for k = 2:num_cell_z+2
            if(abs(dV_X(i,j,k)) < 10^-13)  %USING REAL TAYLOR EXPANSION HERE CAUSES BLOWUP!!!--> WHY?
                Bp_posX(i,j,k) = 1;%1 -  dV_X(i,j,k)/2 + (dV_X(i,j,k))^2/12 - (dV_X(i,j,k))^4/720;
                Bp_negX(i,j,k) =  1;%Bp_posX(i,j)*exp(dV_X(i,j,k));
            end
            if(abs(dV_Y(i,j,k)) < 10^-13)  %using real taylor expansion here works fine...
                Bp_posY(i,j,k) = 1;%1 -  dV_Y(i,j,k)/2 + (dV_Y(i,j,k))^2/12 - (dV_Y(i,j,k))^4/720;
                Bp_negY(i,j,k) = 1;%Bp_posY(i,j,k)*exp(dV_Y(i,j,k));
            end     
            if(abs(dV_Z(i,j,k)) < 10^-13)  %using real taylor expansion here works fine...
                Bp_posZ(i,j,k) = 1;%1 -  dV_Z(i,j,k)/2 + (dV_Z(i,j,k))^2/12 - (dV_Z(i,j,k))^4/720;
                Bp_negZ(i,j,k) = 1;%Bp_posZ(i,j,k)*exp(dV_Z(i,j,k));
            end         
        end
    end
end



%create a structure for this function to return all Bernoulli values
Bernoulli_p_values.Bp_posX = Bp_posX;
Bernoulli_p_values.Bp_negX = Bp_negX;
Bernoulli_p_values.Bp_posY = Bp_posY;
Bernoulli_p_values.Bp_negY = Bp_negY;
Bernoulli_p_values.Bp_posZ = Bp_posZ;
Bernoulli_p_values.Bp_negZ = Bp_negZ;
