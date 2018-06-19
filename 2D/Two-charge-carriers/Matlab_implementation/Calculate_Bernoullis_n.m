function Bernoulli_n_values = Calculate_Bernoullis_n(fullV)
global num_cell 

%This will return a vector/array containing valuess for all 4 Bernoulli p fncs as follows:
%Bernoulli_values = [Bn_posX  Bn_negX Bn_posZ Bn_negZ]

Bn_posX = zeros(num_cell,num_cell);
Bn_negX = zeros(num_cell,num_cell);
Bn_posZ = zeros(num_cell,num_cell);
Bn_negZ = zeros(num_cell,num_cell);
dV_X = zeros(num_cell+1,num_cell+1);
dV_Z = zeros(num_cell+1,num_cell+1);

for i = 2:num_cell+1
    for j = 2:num_cell+1
        dV_X(i,j) = fullV(i,j) - fullV(i-1,j);   %NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
        dV_Z(i,j) = fullV(i,j) - fullV(i,j-1);
    end
end


for i = 1:num_cell
    for j = 1:num_cell
        if(abs(dV_X(i+1,j+1)) < 10^-13)  %to prevent blowup due  to 0 denominator
            Bn_posX(i,j) = 1;% -  dV_X(i+1,j+1)./2. + (dV_X(i+1,j+1))^2./12. - (dV_X(i+1,j+1))^4./720.;
            Bn_negX(i,j) =  Bn_posX(i,j)*exp(dV_X(i+1));  %these will not blow up
        else
            Bn_posX(i,j) = dV_X(i+1,j+1)/(exp(dV_X(i+1,j+1))-1.0);  %the +1's b/c dV's are defined as from 2, and 1 corresponds to the boundary...
            Bn_negX(i,j) =  Bn_posX(i,j)*exp(dV_X(i+1));
        end
        if(abs(dV_Z(i+1,j+1)) < 10^-13)
             Bn_posZ(i,j) = 1;% -  dV_Z(i+1,j+1)./2. + (dV_Z(i+1,j+1))^2./12. - (dV_Z(i+1,j+1))^4./720.;
             Bn_negZ(i,j) =  Bn_posZ(i,j)*exp(dV_Z(i+1,j+1));
        else
            Bn_posZ(i,j) = dV_Z(i+1,j+1)/(exp(dV_Z(i+1,j+1))-1.0);
            Bn_negZ(i,j) =  Bn_posZ(i,j)*exp(dV_Z(i+1,j+1));
        end
    end
end

%create a structure for this function to return all Bernoulli values
Bernoulli_n_values.Bn_posX = Bn_posX;
Bernoulli_n_values.Bn_negX = Bn_negX;
Bernoulli_n_values.Bn_posZ = Bn_posZ;
Bernoulli_n_values.Bn_negZ = Bn_negZ;
