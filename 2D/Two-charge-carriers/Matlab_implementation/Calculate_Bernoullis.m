function Bernoulli_p_values = Calculate_Bernoullis(fullV)
global num_cell l_HTL_int Vt HTL_int_VBstep phi_a l_ETL_int ETL_int_VBstep phi_c BCP_int_VBstep 

%This will return a vector/array containing valuess for all 4 Bernoulli p
%fncs as follows:
%Bernoulli_values = [Bp_posX  Bp_negX Bp_posZ Bp_negZ]

dV_X = zeros(num_cell+1,num_cell+1);
dV_Z = zeros(num_cell+1,num_cell+1);

for i = 2:num_cell+1
    for j = 2:num_cell+1
        dV_X(i,j) = fullV(i,j) - fullV(i-1,j);   %NOTE: all the dV_X are ~0, b/c in x direction I made there be no difference in electric potential
        dV_Z(i,j) = fullV(i,j) - fullV(i,j-1);
    end
end

%ISSUE IS THAT, since my dV_X are near 0 --> in bernoullli fnc get
%denominator: exp^0 - 1 = 0 --> 0 denominator BLOWS UP!!

%NEED TO USE THE FACT THAT THE LIMIT OF x/(e^x -1) as x goeos to 0 is 1
%So if dV is too small --> set bernoulli to = 1...


%add the injection steps here

for i = 1:num_cell
    for j = 1:num_cell
        if(dV_X(i+1,j+1) < 10^-13)  %to prevent blowup due  to 0 denominator
            Bp_posX(i,j) = 1.;
            Bp_negX(i,j) = 1.;
        else
            Bp_posX(i,j) = dV_X(i+1,j+1)/(exp(dV_X(i+1,j+1))-1.0);  %the +1's b/c dV's are defined as from 2, and 1 corresponds to the boundary...
            Bp_negX(i,j) =  Bp_posX(i,j)*exp(dV_X(i+1));
        end
        if(dV_Z(i+1,j+1) < 10^-13)
             Bp_posZ(i,j) = 1.;
             Bp_negZ(i,j) = 1.;
        else
            Bp_posZ(i,j) = dV_Z(i+1,j+1)/(exp(dV_Z(i+1,j+1))-1.0);
            Bp_negZ(i,j) =  Bp_posZ(i,j)*exp(dV_Z(i+1,j+1));
        end
    end
end

%create a structure for this function to return all Bernoulli values
Bernoulli_p_values.Bp_posX = Bp_posX;
Bernoulli_p_values.Bp_negX = Bp_negX;
Bernoulli_p_values.Bp_posZ = Bp_posZ;
Bernoulli_p_values.Bp_negZ = Bp_negZ;