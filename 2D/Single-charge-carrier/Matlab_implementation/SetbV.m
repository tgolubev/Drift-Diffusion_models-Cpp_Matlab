function bV = SetbV(p_matrix, epsilon)
global V_leftBC V_bottomBC V_topBC V_rightBC N CV;

 %set up rhs of Poisson equation. Note for epsilons, are assuming that
    %epsilons at the boundaries are the same as espilon 1 cell into interior of
    %device
    index = 0;
    for j = 1:N
        if(j ==1)  %different for 1st subblock
            for i = 1:N
                index = index +1;
                if (i==1)  %1st element has 2 BC's
                    bV(index,1) = CV*(p_matrix(i,j)) + epsilon(i,j)*(V_leftBC(1) + V_bottomBC);   %RECALL matrix has +2 down diagonals, sign flipped from 1D version
                elseif (i==N)
                    bV(index,1) = CV*(p_matrix(i,j)) + epsilon(i,j)*(V_rightBC(1) + V_bottomBC);
                else
                    bV(index,1) = CV*(p_matrix(i,j)) + epsilon(i,j)*V_bottomBC;
                end
            end
        elseif(j == N)  %different for last subblock
            for i = 1:N
                index = index +1;
                if (i==1)  %1st element has 2 BC's
                    bV(index,1) = CV*(p_matrix(i,j)) + epsilon(i,j)*(V_leftBC(N) + V_topBC);
                elseif (i==N)
                    bV(index,1) = CV*(p_matrix(i,j)) + epsilon(i,j)*(V_rightBC(N) + V_topBC);
                else
                    bV(index,1) = CV*(p_matrix(i,j)) + epsilon(i,j)*V_topBC;
                end
            end
        else %interior subblocks
            for i = 1:N
                index = index +1;
                if(i==1)
                    bV(index,1) = CV*(p_matrix(i,j)) + epsilon(i,j)*V_leftBC(j);
                elseif(i==N)
                    bV(index,1) = CV*(p_matrix(i,j)) + epsilon(i,j)*V_rightBC(j);
                else
                    bV(index,1) = CV*(p_matrix(i,j));
                end
            end
        end
    end
