function bp = Setbp_2D(Bernoulli_p_values, p_mob, Up)
global Cp num_elements p_topBC p_leftBC p_rightBC p_bottomBC N;

% bp = zeros(num_elements,1);
%Note: assuming for BC's that mobility  on bndry is same as mobility just
%inside the boundary (so not doing any mobs averaging there)

%extract variables from struct (for brevity in eqns)
Bp_posX = Bernoulli_p_values.Bp_posX;
Bp_negX = Bernoulli_p_values.Bp_negX;
Bp_posZ = Bernoulli_p_values.Bp_posZ;
Bp_negZ = Bernoulli_p_values.Bp_negZ;

%calculate main part here
bp = Cp*Up;

%add on BC's
index = 0;
    for j = 1:N
        if(j ==1)  %different for 1st sub-block
            for i = 1:N   
                index = index +1;
                if (i==1)  %1st element has 2 BC's  %NOTE: p_mob are +1, b/c corresopnd to FULL device while Up and bp's correspond to interior elements only
                    %Bp are +1, because defined from 2....
                    bp(index,1) = bp(index,1) + p_mob(i+1,j+1)*(Bp_posX(i+1,j+1)*p_leftBC(1) + Bp_posZ(i+1,j+1)*p_bottomBC);  %NOTE: rhs is +Cp*Up, b/c diagonal elements are + here, flipped sign from 1D version              
                elseif (i==N)
                    bp(index,1) = bp(index,1) + p_mob(i+1,j+1)*(Bp_posZ(i+1,j+1)*p_bottomBC + Bp_negX(i+1+1,j+1)*p_rightBC(1));
                else
                    bp(index,1) = bp(index,1) + p_mob(i+1,j+1)*Bp_posZ(i+1,j+1)*p_bottomBC;
                end              
            end
        elseif(j == N)  %different for last subblock
            for i = 1:N   
                index = index +1;
                if (i==1)  %1st element has 2 BC's
                    bp(index,1) = bp(index,1) + p_mob(i+1,j+1)*(Bp_posX(i+1,j+1)*p_leftBC(N) + Bp_negZ(i+1,j+1+1)*p_topBC);
                elseif (i==N)
                    bp(index,1) = bp(index,1) + p_mob(i+1,j+1)*(Bp_negX(i+1+1,j+1)*p_rightBC(N) + Bp_negZ(i+1,j+1+1)*p_topBC);
                else
                    bp(index,1) = bp(index,1) + p_mob(i+1,j+1)*Bp_negZ(i+1,j+1+1)*p_topBC;
                end
            end
        else %interior subblocks
            for i = 1:N  
                index = index +1;
                if(i==1)
                    bp(index,1) = bp(index,1) + p_mob(i+1,j+1)*Bp_posX(i+1,j+1)*p_leftBC(j);
                elseif(i==N)
                    bp(index,1) = bp(index,1) + p_mob(i+1,j+1)*Bp_negX(i+1+1,j+1)*p_rightBC(j);
                end
            end
        end
    end
