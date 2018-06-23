function Bernoulli_p_values = Calculate_Bernoullis(fullV)
global num_cell 

%This will return a vector/array containing valuess for all 4 Bernoulli p fncs as follows:
%Bernoulli_values = [Bp_posX  Bp_negX Bp_posZ Bp_negZ]

Bp_posX = zeros(num_cell+1,num_cell+1);
Bp_negX = zeros(num_cell+1,num_cell+1);
Bp_posZ = zeros(num_cell+1,num_cell+1);
Bp_negZ = zeros(num_cell+1,num_cell+1);
dV_X = zeros(num_cell+1,num_cell+1);
dV_Z = zeros(num_cell+1,num_cell+1);

for i = 2:num_cell+1
    for j = 2:num_cell+1
        dV_X(i,j) = fullV(i,j) - fullV(i-1,j);   %NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
        dV_Z(i,j) = fullV(i,j) - fullV(i,j-1);
    end
end


for i = 2:num_cell+1
    for j = 2:num_cell+1
        if(abs(dV_X(i,j)) < 10^-13)  %to prevent blowup due  to 0 denominator
            %taylor expand x/(e^x -1)
            Bp_posX(i,j) = 1;%1 -  dV_X(i+1,j+1)/2 + (dV_X(i+1,j+1))^2/12 - (dV_X(i+1,j+1))^4/720;
            Bp_negX(i,j) =  1;%Bp_posX(i,j)*exp(dV_X(i+1,j+1));
        else
            Bp_posX(i,j) = 1;%dV_X(i+1,j+1)/(exp(dV_X(i+1,j+1))-1.0);  %the +1's b/c dV's are defined as from 2, and 1 corresponds to the boundary...
            Bp_negX(i,j) = 1;%Bp_posX(i,j)*exp(dV_X(i,j));
        end
        if(abs(dV_Z(i,j)) < 10^-13)
              Bp_posZ(i,j) = 1;%1 -  dV_X(i+1,j+1)/2 + (dV_X(i+1,j+1))^2/12 - (dV_X(i+1,j+1))^4/720;
              Bp_negZ(i,j) =  1;%Bp_posZ(i,j)*exp(dV_Z(i+1,j+1));
        else
            Bp_posZ(i,j) = dV_Z(i,j)/(exp(dV_Z(i,j))-1.0);
            Bp_negZ(i,j) =  Bp_posZ(i,j)*exp(dV_Z(i,j));
         end
    end
end

%create a structure for this function to return all Bernoulli values
Bernoulli_p_values.Bp_posX = Bp_posX;
Bernoulli_p_values.Bp_negX = Bp_negX;
Bernoulli_p_values.Bp_posZ = Bp_posZ;
Bernoulli_p_values.Bp_negZ = Bp_negZ;
