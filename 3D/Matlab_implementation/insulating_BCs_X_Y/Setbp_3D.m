function bp = Setbp_3D(Bernoulli_p_values, p_mob, Up)
global Cp num_elements p_topBC p_leftBC_x p_rightBC_x p_leftBC_y p_rightBC_y p_bottomBC N;

% bp = zeros(num_elements,1);
%Note: assuming for BC's that mobility  on bndry is same as mobility just
%inside the boundary (so not doing any mobs averaging there)

%extract variables from struct (for brevity in eqns)
Bp_posX = Bernoulli_p_values.Bp_posX;
Bp_negX = Bernoulli_p_values.Bp_negX;
Bp_posY = Bernoulli_p_values.Bp_posY;
Bp_negY = Bernoulli_p_values.Bp_negY;
Bp_posZ = Bernoulli_p_values.Bp_posZ;
Bp_negZ = Bernoulli_p_values.Bp_negZ;

%calculate main part here
bp = Cp*Up;

%add on BC's
index = 0;
for k = 1:N
    if (k == 1)  %--------------------------------------------------------------------
        for j = 1:N
            if(j == 1)  %different for 1st subblock
                for i = 1:N
                    index = index +1;
                    if (i == 1)  %1st element has 2 BC's
                        bp(index,1) = bp(index,1) + p_mob(i+1,j+1,k+1)*p_leftBC_x(1,1)*Bp_posX(i+1,j+1,k+1) + p_mob(i+1,j+1,k+1)*p_leftBC_y(1,1)*Bp_posY(i+1,j+1,k+1) + p_mob(i+1,j+1,k+1)*p_bottomBC*Bp_posZ(i+1,j+1,k+1);  %+1 b/c netcharge and p_mob include endpoints but i,j index only the interior
                    elseif (i==N)
                        bp(index,1) = bp(index,1) + p_mob(N+1+1, 1+1, 1+1)*p_rightBC_x(1,1)*Bp_negX(i+1+1,j+1,k+1) + p_mob(N+1, 0+1, 1+1)*p_leftBC_y(N,1)*Bp_posY(i+1,j+1,k+1) + p_mob(N+1, 1+1, 0+1)*p_bottomBC*Bp_posZ(i+1,j+1,k+1);
                    else  %middle elements in 1st subblock
                        bp(index,1) = bp(index,1) + p_mob(i+1, 0+1, 1+1)*p_leftBC_y(i, 1)*Bp_posY(i+1,j+1,k+1) + p_mob(i+1, 1+1, 0+1)*p_bottomBC*Bp_posZ(i+1,j+1,k+1);
                    end
                end
            elseif(j == N)  %different for last subblock within k=1 subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)  %1st element has 2 BC's
                        bp(index,1) = bp(index,1) + p_mob(0+1,N+1, 1+1)*p_leftBC_x(N,1)*Bp_posX(i+1,j+1,k+1) + p_mob(1+1, N+1+1, 1+1)*p_rightBC_y(1,1)*Bp_negY(i+1,j+1+1,k+1) + p_mob(1+1, N+1, 0+1)*p_bottomBC*Bp_posZ(i+1,j+1,k+1);
                    elseif (i==N)
                        bp(index,1) = bp(index,1) + p_mob(N+1+1,N+1, 1+1)*p_rightBC_x(N,1)*Bp_negX(i+1+1,j+1,k+1) + p_mob(N+1, N+1+1, 1+1)*p_rightBC_y(N,1)*Bp_negY(i+1,j+1+1,k+1) + p_mob(N+1,N+1, 0+1)*p_bottomBC*Bp_posZ(i+1,j+1,k+1);
                    else %inner rows of Nth subblock
                        bp(index,1) = bp(index,1) + p_mob(i+1,N+1+1, 1+1)*p_rightBC_y(i,1)*Bp_negY(i+1,j+1+1,k+1) + p_mob(i+1, N+1, 0+1)*p_bottomBC*Bp_posZ(i+1,j+1,k+1);
                    end
                end
            else  %interior subblocks of 1st (k = 1) subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)
                        bp(index,1) = bp(index,1) + p_mob(0+1, j+1, 1+1)*p_leftBC_x(j,1)*Bp_posX(i+1,j+1,k+1) + p_mob(1+1, j+1, 0+1)*p_bottomBC*Bp_posZ(i+1,j+1,k+1);
                    elseif (i==N)
                        bp(index,1) = bp(index,1) + p_mob(N+1+1, j+1, 1+1)*p_rightBC_x(j,1)*Bp_negX(i+1+1,j+1,k+1) + p_mob(N+1, j+1, 0+1)*p_bottomBC*Bp_posZ(i+1,j+1,k+1);
                    else
                        bp(index,1) = bp(index,1) + p_mob(i+1, j+1, 0+1)*p_bottomBC*Bp_posZ(i+1,j+1,k+1);
                    end
                end
            end
        end
    elseif (k == N)  %last subblock group
        for j = 1:N
            if(j == 1)  %different for 1st subblock
                for i = 1:N
                    index = index +1;
                    if (i == 1)  %1st element has 2 BC's
                        bp(index,1) = bp(index,1) + p_mob(0+1, 1+1, N+1)*p_leftBC_x(1,N)*Bp_posX(i+1,j+1,k+1) + p_mob(1+1, 0+1, N+1)*p_leftBC_y(1,N)*Bp_posY(i+1,j+1,k+1) + p_mob(1+1, 1+1, N+1+1)*p_topBC*Bp_negZ(i+1,j+1,k+1+1);  %+1 b/c netcharge and p_mob include endpoints but i,j index only the interior
                    elseif (i==N)
                        bp(index,1) = bp(index,1) + p_mob(N+1+1, 1+1, N+1)*p_rightBC_x(1,N)*Bp_negX(i+1+1,j+1,k+1) + p_mob(N+1, 0+1, N+1)*p_leftBC_y(N,N)*Bp_posY(i+1,j+1,k+1) + p_mob(N+1, 1+1, N+1+1)*p_topBC*Bp_negZ(i+1,j+1,k+1+1);
                    else  %middle elements in 1st subblock
                        bp(index,1) = bp(index,1) + p_mob(i+1, 0+1, N+1)*p_leftBC_y(i, N)*Bp_posY(i+1,j+1,k+1) + p_mob(i+1, 1+1, N+1+1)*p_topBC*Bp_negZ(i+1,j+1,k+1+1);
                    end
                end
            elseif(j == N)  %different for last subblock within k=N subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)  %1st element has 2 BC's
                        bp(index,1) = bp(index,1) + p_mob(0+1,N+1, N+1)*p_leftBC_x(N,N)*Bp_posX(i+1,j+1,k+1) + p_mob(1+1, N+1+1, N+1)*p_rightBC_y(1,N)*Bp_negY(i+1,j+1+1,k+1) + p_mob(1+1, N+1, N+1+1)*p_topBC*Bp_negZ(i+1,j+1,k+1+1);
                    elseif (i==N)
                        bp(index,1) = bp(index,1) + p_mob(N+1+1,N+1, N+1)*p_rightBC_x(N,N)*Bp_negX(i+1+1,j+1,k+1) + p_mob(N+1, N+1+1, N+1)*p_rightBC_y(N,N)*Bp_negY(i+1,j+1+1,k+1) + p_mob(N+1,N+1, N+1+1)*p_topBC*Bp_negZ(i+1,j+1,k+1+1);
                    else %inner rows of Nth subblock
                        bp(index,1) = bp(index,1) + p_mob(i+1,N+1+1, N+1)*p_rightBC_y(i,N)*Bp_negY(i+1,j+1+1,k+1) + p_mob(i+1, N+1, N+1+1)*p_topBC*Bp_negZ(i+1,j+1,k+1+1);
                    end
                end
            else  %interior subblocks of last (k = N) subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)
                        bp(index,1) = bp(index,1) + p_mob(0+1, j+1, N+1)*p_leftBC_x(j,N)*Bp_posX(i+1,j+1,k+1) + p_mob(1+1, j+1, N+1+1)*p_topBC*Bp_negZ(i+1,j+1,k+1+1);
                    elseif (i==N)
                        bp(index,1) = bp(index,1) + p_mob(N+1+1, j+1, N+1)*p_rightBC_x(j,N)*Bp_negX(i+1+1,j+1,k+1) + p_mob(N+1, j+1, N+1+1)*p_topBC*Bp_negZ(i+1,j+1,k+1+1);
                    else
                        bp(index,1) = bp(index,1) + p_mob(i+1, j+1, N+1+1)*p_topBC*Bp_negZ(i+1,j+1,k+1+1);
                    end
                end
            end
        end
    else   %interior subblock groups (k=2:N-1)
        for j = 1:N
            if(j == 1)  %different for 1st subblock
                for i = 1:N
                    index = index +1;
                    if (i == 1)  %1st element has 2 BC's
                        bp(index,1) = bp(index,1) + p_mob(0+1, 1+1, k+1)*p_leftBC_x(1,k)*Bp_posX(i+1,j+1,k+1)  + p_mob(1+1, 0+1, k+1)*p_leftBC_y(1,k)*Bp_posY(i+1,j+1,k+1);  %+1 b/c netcharge and p_mob include endpoints but i,j index only the interior
                    elseif (i==N)
                        bp(index,1) = bp(index,1) + p_mob(N+1+1, 1+1, k+1)*p_rightBC_x(1,k)*Bp_negX(i+1+1,j+1,k+1) + p_mob(N+1, 0+1, k+1)*p_leftBC_y(N,k)*Bp_posY(i+1,j+1,k+1);
                    else  %middle elements in 1st subblock
                        bp(index,1) = bp(index,1) + p_mob(i+1, 0+1, k+1)*p_leftBC_y(i, k)*Bp_posY(i+1,j+1,k+1);
                    end
                end
            elseif(j == N)  %different for last subblock
                for i = 1:N
                    index = index +1;
                    if (i==1)  %1st element has 2 BC's
                        bp(index,1) = bp(index,1) + p_mob(0+1,N+1, k+1)*p_leftBC_x(N,k)*Bp_posX(i+1,j+1,k+1) + p_mob(1+1, N+1+1, k+1)*p_rightBC_y(1,k)*Bp_negY(i+1,j+1+1,k+1);
                    elseif (i==N)
                        bp(index,1) = bp(index,1) + p_mob(N+1+1,N+1, k+1)*p_rightBC_x(N,k)*Bp_negX(i+1+1,j+1,k+1) + p_mob(N+1, N+1+1, k+1)*p_rightBC_y(N,k)*Bp_negY(i+1,j+1+1,k+1);
                    else %inner rows of Nth subblock
                        bp(index,1) = bp(index,1) + p_mob(i+1,N+1+1, k+1)*p_rightBC_y(i,k)*Bp_negY(i+1,j+1+1,k+1);
                    end
                end
            else  %interior subblocks
                for i = 1:N
                    index = index +1;
                    if (i==1)
                        bp(index,1) = bp(index,1) + p_mob(0+1, j+1, k+1)*p_leftBC_x(j,k)*Bp_posX(i+1,j+1,k+1) ;
                    elseif (i==N)
                        bp(index,1) = bp(index,1) + p_mob(N+1+1, j+1, k+1)*p_rightBC_x(j,k)*Bp_negX(i+1+1,j+1,k+1);
                    end
                end
            end
        end
    end
end
