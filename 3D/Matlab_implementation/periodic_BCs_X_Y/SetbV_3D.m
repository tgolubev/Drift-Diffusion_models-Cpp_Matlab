function bV = SetbV_3D(p, epsilon)
global V_bottomBC V_topBC N CV;

 %set up rhs of Poisson equation. Note for epsilons, are assuming that
    %epsilons at the boundaries are the same as espilon cell into interior of
    %device
    
bV = CV*(p);  %note: this are in vector form

index = 0;
for i = 1:N+1
        for j = 1:N+1
            for k = 1:N
                index = index + 1;
                if (k == 1)  %bottom BC                  
                    bV(index) = bV(index) + epsilon(i+1,j+1,0+1)*V_bottomBC(i,j);
                end
                %middle elements have no BC's
                if (k == N) %top BC
                    bV(index) = bV(index) + epsilon(i+1,j+1,N+1)*V_topBC(i,j);
                end
            end
        end
end
                        
                   
%Notice how much simpler rhs is with this ordering
%and including PBC's in the matrix!