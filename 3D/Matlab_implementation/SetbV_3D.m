function bV = SetbV_3D(p_matrix, n_matrix, epsilon)
global V_leftBC_x V_leftBC_y V_bottomBC V_topBC V_rightBC_x V_rightBC_y N CV;

 %set up rhs of Poisson equation. Note for epsilons, are assuming that
    %epsilons at the boundaries are the same as espilon cell into interior of
    %device
    
   
bV = CV*(p_matrix-n_matrix);

index = 0;
for k = 1:N
    if (k == 1)  %--------------------------------------------------------------------
        for j = 1:N
            if(j == 1)  %different for 1st subblock
                for i = 1:N
                    index = index +1;
                    if (i == 1)  %1st element has 2 BC's
                        bV(index,1) = bV(index,1) + epsilon(0+1, 1+1, 1+1)*V_leftBC_x(1,1) + epsilon(1+1, 0+1, 1+1)*V_leftBC_y(1,1)  + epsilon(1+1, 1+1, 0+1)*V_bottomBC;  %+1 b/c netcharge and epsilon include endpoints but i,j index only the interior
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, 1+1, 1+1)*V_rightBC_x(1,1) + epsilon(N+1, 0+1, 1+1)*V_leftBC_y(N,1) + epsilon(N+1, 1+1, 0+1)*V_bottomBC;
                    else  %middle elements in 1st subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1, 0+1, 1+1)*V_leftBC_y(i, 1) + epsilon(i+1, 1+1, 0+1)*V_bottomBC;
                    end
                end
            elseif(j == N)  %different for last subblock within k=1 subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)  %1st element has 2 BC's
                        bV(index,1) = bV(index,1) + epsilon(0+1,N+1, 1+1)*V_leftBC_x(N,1) + epsilon(1+1, N+1+1, 1+1)*V_rightBC_y(1,1) + epsilon(1+1, N+1, 0+1)*V_bottomBC;
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1,N+1, 1+1)*V_rightBC_x(N,1) + epsilon(N+1, N+1+1, 1+1)*V_rightBC_y(N,1) + epsilon(N+1,N+1, 0+1)*V_bottomBC;
                    else %inner rows of Nth subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1,N+1+1, 1+1)*V_rightBC_y(i,1) + epsilon(i+1, N+1, 0+1)*V_bottomBC;
                    end
                end
            else  %interior subblocks of 1st (k = 1) subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)
                        bV(index,1) = bV(index,1) + epsilon(0+1, j+1, 1+1)*V_leftBC_x(j,1) + epsilon(1+1, j+1, 0+1)*V_bottomBC;
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, j+1, 1+1)*V_rightBC_x(j,1) + epsilon(N+1, j+1, 0+1)*V_bottomBC;
                    else
                        bV(index,1) = bV(index,1) + epsilon(i+1, j+1, 0+1)*V_bottomBC;
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
                        bV(index,1) = bV(index,1) + epsilon(0+1, 1+1, N+1)*V_leftBC_x(1,N) + epsilon(1+1, 0+1, N+1)*V_leftBC_y(1,N)  + epsilon(1+1, 1+1, N+1+1)*V_topBC;  %+1 b/c netcharge and epsilon include endpoints but i,j index only the interior
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, 1+1, N+1)*V_rightBC_x(1,N) + epsilon(N+1, 0+1, N+1)*V_leftBC_y(N,N) + epsilon(N+1, 1+1, N+1+1)*V_topBC;
                    else  %middle elements in 1st subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1, 0+1, N+1)*V_leftBC_y(i, N) + epsilon(i+1, 1+1, N+1+1)*V_topBC;
                    end
                end
            elseif(j == N)  %different for last subblock within k=N subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)  %1st element has 2 BC's
                        bV(index,1) = bV(index,1) + epsilon(0+1,N+1, N+1)*V_leftBC_x(N,N) + epsilon(1+1, N+1+1, N+1)*V_rightBC_y(1,N) + epsilon(1+1, N+1, N+1+1)*V_topBC;
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1,N+1, N+1)*V_rightBC_x(N,N) + epsilon(N+1, N+1+1, N+1)*V_rightBC_y(N,N) + epsilon(N+1,N+1, N+1+1)*V_topBC;
                    else %inner rows of Nth subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1,N+1+1, N+1)*V_rightBC_y(i,N) + epsilon(i+1, N+1, N+1+1)*V_topBC;
                    end
                end
            else  %interior subblocks of last (k = N) subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)
                        bV(index,1) = bV(index,1) + epsilon(0+1, j+1, N+1)*V_leftBC_x(j,N) + epsilon(1+1, j+1, N+1+1)*V_topBC;
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, j+1, N+1)*V_rightBC_x(j,N) + epsilon(N+1, j+1, N+1+1)*V_topBC;
                    else
                        bV(index,1) = bV(index,1) + epsilon(i+1, j+1, N+1+1)*V_topBC;
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
                        bV(index,1) = bV(index,1) + epsilon(0+1, 1+1, k+1)*V_leftBC_x(1,k) + epsilon(1+1, 0+1, k+1)*V_leftBC_y(1,k);  %+1 b/c netcharge and epsilon include endpoints but i,j index only the interior
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, 1+1, k+1)*V_rightBC_x(1,k) + epsilon(N+1, 0+1, k+1)*V_leftBC_y(N,k);
                    else  %middle elements in 1st subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1, 0+1, k+1)*V_leftBC_y(i, k);
                    end
                end
            elseif(j == N)  %different for last subblock
                for i = 1:N
                    index = index +1;
                    if (i==1)  %1st element has 2 BC's
                        bV(index,1) = bV(index,1) + epsilon(0+1,N+1, k+1)*V_leftBC_x(N,k) + epsilon(1+1, N+1+1, k+1)*V_rightBC_y(1,k);
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1,N+1, k+1)*V_rightBC_x(N,k) + epsilon(N+1, N+1+1, k+1)*V_rightBC_y(N,k);
                    else %inner rows of Nth subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1,N+1+1, k+1)*V_rightBC_y(i,k);
                    end
                end
            else  %interior subblocks
                for i = 1:N
                    index = index +1;
                    if (i==1)
                        bV(index,1) = bV(index,1) + epsilon(0+1, j+1, k+1)*V_leftBC_x(j,k);
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, j+1, k+1)*V_rightBC_x(j,k);
                    end
                end
            end
        end
    end
end
