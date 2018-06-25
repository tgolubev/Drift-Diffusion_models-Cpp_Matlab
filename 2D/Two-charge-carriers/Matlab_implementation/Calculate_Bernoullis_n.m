function Bernoulli_n_values = Calculate_Bernoullis_n(fullV)
global num_cell 

%This will return a vector/array containing valuess for all 4 Bernoulli p fncs as follows:
%Bernoulli_values = [Bn_posX  Bn_negX Bn_posZ Bn_negZ]
%matlab suggests to not preallocate these--> b/c are output of the
%function. Already allocated when call the function.

% Bn_posX = zeros(num_cell+1,num_cell+1);
% Bn_negX = zeros(num_cell+1,num_cell+1);
% Bn_posZ = zeros(num_cell+1,num_cell+1);
% Bn_negZ = zeros(num_cell+1,num_cell+1);
dV_X = zeros(num_cell+1,num_cell+1);
dV_Z = zeros(num_cell+1,num_cell+1);

for i = 2:num_cell+1
    for j = 2:num_cell+1
        dV_X(i,j) = fullV(i,j) - fullV(i-1,j);   %NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
        dV_Z(i,j) = fullV(i,j) - fullV(i,j-1);
    end
end

% matrix operations for Bernoulli's, SAVES SIGNIFICANT CPU then overwrite the elements 
%which need overwriting...
%BE CAREFUL!, NEED TO USE THE . to indicate element-by-element operations
Bn_posX = dV_X./(exp(dV_X) - 1.0);
Bn_negX = Bn_posX.*exp(dV_X);

Bn_posZ = dV_Z./(exp(dV_Z) - 1.0);
Bn_negZ = Bn_posZ.*exp(dV_Z);


for i = 2:num_cell+1
    for j = 2:num_cell+1
        if(abs(dV_X(i,j)) < 10^-13)  %to prevent blowup due  to 0 denominator
            Bn_posX(i,j) = 1;%1 -  dV_X(i+1,j+1)/2 + (dV_X(i+1,j+1))^2/12 - (dV_X(i+1,j+1))^4/720;
            Bn_negX(i,j) =  1;%Bn_posX(i,j)*exp(dV_X(i+1));  %these will not blow up
        end
%         else
%             Bn_posX(i,j) = 1;%dV_X(i,j)/(exp(dV_X(i,j))-1.0);  %the +1's b/c dV's are defined as from 2, and 1 corresponds to the boundary...
%             Bn_negX(i,j) =  1;%Bn_posX(i,j)*exp(dV_X(i,j));
%         end
        if(abs(dV_Z(i,j)) < 10^-13)
             Bn_posZ(i,j) = 1;%1 -  dV_X(i+1,j+1)/2 + (dV_X(i+1,j+1))^2/12 - (dV_X(i+1,j+1))^4/720;
             Bn_negZ(i,j) =  1;%Bn_posZ(i,j)*exp(dV_Z(i+1,j+1));
        end
%         else
%             Bn_posZ(i,j) = dV_Z(i,j)/(exp(dV_Z(i,j))-1.0);
%             Bn_negZ(i,j) =  Bn_posZ(i,j)*exp(dV_Z(i,j));
%         end
    end
end

%create a structure for this function to return all Bernoulli values
Bernoulli_n_values.Bn_posX = Bn_posX;
Bernoulli_n_values.Bn_negX = Bn_negX;
Bernoulli_n_values.Bn_posZ = Bn_posZ;
Bernoulli_n_values.Bn_negZ = Bn_negZ;
