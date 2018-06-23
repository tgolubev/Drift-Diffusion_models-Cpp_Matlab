function bn = Setbn_2D(Bernoulli_n_values, n_mob, Un)
global Cn num_elements n_topBC n_leftBC n_rightBC n_bottomBC N;

bn = zeros(num_elements,1);

%extract variables from struct (for brevity in eqns)
Bn_posX = Bernoulli_n_values.Bn_posX;
Bn_negX = Bernoulli_n_values.Bn_negX;
Bn_posZ = Bernoulli_n_values.Bn_posZ;
Bn_negZ = Bernoulli_n_values.Bn_negZ;

index = 0;
    for j = 1:N
        if(j ==1)  %different for 1st subblock
            for i = 1:N
                index = index +1;
                if (i==1)  %1st element has 2 BC's
                    bn(index,1) = Cn*Un(i,j) + n_mob(i+1,j+1)*(Bn_negX(i+1,j+1)*n_leftBC(1) + Bn_negZ(i+1,j+1)*n_bottomBC);  %NOTE: rhs is +Cp*Un, b/c diagonal elements are + here, flipped sign from 1D version                   
                elseif (i==N)
                    bn(index,1) = Cn*Un(i,j) + n_mob(i+1,j+1)*(Bn_negZ(i+1,j+1)*n_bottomBC + Bn_posX(i+1+1,j+1)*n_rightBC(1));
                else
                    bn(index,1) = Cn*Un(i,j) + n_mob(i,j)*Bn_negZ(i+1,j+1)*n_bottomBC;
                end
            end
        elseif(j == N)  %different for last subblock
            for i = 1:N
                index = index +1;
                if (i==1)  %1st element has 2 BC's
                    bn(index,1) = Cn*Un(i,j) + n_mob(i+1,j+1)*(Bn_negX(i+1,j+1)*n_leftBC(N) + Bn_posZ(i+1,j+1+1)*n_topBC);
                elseif (i==N)
                    bn(index,1) = Cn*Un(i,j) + n_mob(i+1,j+1)*(Bn_posX(i+1+1,j+1)*n_rightBC(N) + Bn_posZ(i+1,j+1+1)*n_topBC);
                else
                    bn(index,1) = Cn*Un(i,j) + n_mob(i+1,j+1)*Bn_posZ(i+1,j+1+1)*n_topBC;
                end
            end
        else %interior subblocks
            for i = 1:N
                index = index +1;
                if(i==1)
                    bn(index,1) = Cn*Un(i,j) + n_mob(i+1,j+1)*Bn_negX(i+1,j+1)*n_leftBC(j);
                elseif(i==N)
                    bn(index,1) = Cn*Un(i,j) + n_mob(i+1,j+1)*Bn_posX(i+1+1,j+1)*n_rightBC(j);
                else
                    bn(index,1) = Cn*Un(i,j);
                end
            end
        end
    end