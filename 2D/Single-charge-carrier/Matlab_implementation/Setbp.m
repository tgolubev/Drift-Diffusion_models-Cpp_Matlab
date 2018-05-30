function bp = Setbp(Bernoulli_p_values, p_mob, Up)
global Cp num_elements p_topBC p_leftBC p_rightBC p_bottomBC N;

bp = zeros(num_elements,1);
%Note: assuming for BC's that mobility  on bndry is same as mobility just
%inside the boundary (so not doing any mobs averaging)

%extract variables from struct (for brevity in eqns)
Bp_posX = Bernoulli_p_values.Bp_posX;
Bp_negX = Bernoulli_p_values.Bp_negX;
Bp_posZ = Bernoulli_p_values.Bp_posZ;
Bp_negZ = Bernoulli_p_values.Bp_negZ;

index = 0;
    for j = 1:N
        if(j ==1)  %different for 1st subblock
            for i = 1:N   
                index = index +1;
                if (i==1)  %1st element has 2 BC's
                    bp(index,1) = Cp*Up(i,j) + p_mob(i,j)*(Bp_posX(i,j)*p_leftBC(1) + Bp_posZ(i,j)*p_bottomBC);  %NOTE: rhs is +Cp*Up, b/c diagonal elements are + here, flipped sign from 1D version              
                elseif (i==N)
                    bp(index,1) = Cp*Up(i,j) + p_mob(i,j)*(Bp_posZ(i,j)*p_bottomBC + Bp_negX(i+1,j)*p_rightBC(1));
                else
                    bp(index,1) = Cp*Up(i,j) + p_mob(i,j)*Bp_posZ(i,j)*p_bottomBC;
                end              
            end
        elseif(j == N)  %different for last subblock
            for i = 1:N   
                index = index +1;
                if (i==1)  %1st element has 2 BC's
                    bp(index,1) = Cp*Up(i,j) + p_mob(i,j)*(Bp_posX(i,j)*p_leftBC(N) + Bp_negZ(i,j+1)*p_topBC);
                elseif (i==N)
                    bp(index,1) = Cp*Up(i,j) + p_mob(i,j)*(Bp_negX(i+1,j)*p_rightBC(N) + Bp_negZ(i,j+1)*p_topBC);
                else
                    bp(index,1) = Cp*Up(i,j) + p_mob(i,j)*Bp_negZ(i,j+1)*p_topBC;
                end
            end
        else %interior subblocks
            for i = 1:N  
                index = index +1;
                if(i==1)
                    bp(index,1) = Cp*Up(i,j) + p_mob(i,j)*Bp_posX(i,j)*p_leftBC(j);
                elseif(i==N)
                    bp(index,1) = Cp*Up(i,j) + p_mob(i,j)*Bp_negX(i+1,j)*p_rightBC(j);
                else
                    bp(index,1) = Cp*Up(i,j);
                end
            end
        end
    end