function Bernoulli_p_values = Calculate_Bernoullis_p(fullV)
global num_cell 

%This will return a vector/array containing valuess for all 4 Bernoulli p fncs as follows:
%Bernoulli_values = [Bp_posX  Bp_negX Bp_posZ Bp_negZ]

%matlab suggests to not preallocate these--> b/c are output of the
%function. Already allocated when call the function.
% Bp_posX = zeros(num_cell+1,num_cell+1);
% Bp_negX = zeros(num_cell+1,num_cell+1);
% Bp_posZ = zeros(num_cell+1,num_cell+1);
% Bp_negZ = zeros(num_cell+1,num_cell+1);
dV_X = zeros(num_cell+1,num_cell+1,num_cell+1);
dV_Y = zeros(num_cell+1,num_cell+1,num_cell+1);
dV_Z = zeros(num_cell+1,num_cell+1,num_cell+1);

for i = 2:num_cell+1
    for j = 2:num_cell+1
        for k = 2:num_cell+1
            dV_X(i,j,k) = fullV(i,j,k) - fullV(i-1,j,k);   %NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
            dV_Y(i,j,k) = fullV(i,j,k) - fullV(i,j-1,k);
            dV_Z(i,j,k) = fullV(i,j,k) - fullV(i,j,k-1);
        end
    end
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

for i = 2:num_cell+1
    for j = 2:num_cell+1
        for k = 2:num_cell+1
            if(abs(dV_X(i,j,k)) < 10^-13)  %to prevent blowup due  to 0 denominator
                %taylor expand x/(e^x -1)
                Bp_posX(i,j,k) = 1;%1 -  dV_X(i,j)/2 + (dV_X(i,j))^2/12 - (dV_X(i+1,j+1))^4/720;
                Bp_negX(i,j,k) =  1;%Bp_posX(i,j)*exp(dV_X(i+1,j+1));
            end
            if(abs(dV_Y(i,j,k)) < 10^-13)
                Bp_posY(i,j,k) = 1;%1 -  dV_X(i,j)/2 + (dV_X(i,j))^2/12 - (dV_X(i,j))^4/720;
                Bp_negY(i,j,k) =  1;%Bp_posZ(i,j)*exp(dV_Z(i,j));
            end     
            if(abs(dV_Z(i,j,k)) < 10^-13)
                Bp_posZ(i,j,k) = 1;%1 -  dV_X(i,j)/2 + (dV_X(i,j))^2/12 - (dV_X(i,j))^4/720;
                Bp_negZ(i,j,k) =  1;%Bp_posZ(i,j)*exp(dV_Z(i,j));
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
