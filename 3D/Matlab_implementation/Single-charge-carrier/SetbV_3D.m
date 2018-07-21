function bV = SetbV_3D(p, epsilon)
global V_bottomBC V_topBC Nx Ny Nz CV;

 %set up rhs of Poisson equation. Note for epsilons, are assuming that
    %epsilons at the boundaries are the same as espilon cell into interior of
    %device 
    
bV = CV*(p)./18.; %try scaling....  %note: this are in vector form

index = 0;
for i = 1:Nx+1
    for j = 1:Ny+1
        
        for k = 1:Nz+1
            index = index + 1;
            if (k == 1)  %bottom BC
                bV(index) = bV(index) + (epsilon(i+1,j+1,0+1)./18.)*V_bottomBC(i,j);
            end
            %middle elements have no BC's
            if (k == Nz+1) %top BC\\
                
                %1st change the charge part everywhere
                
%                 bV(index) = bV(index)/2;   %this is due to I used -epsilon_Z instead of -2*epsilon_Z in Poisson for neuman BC
                %NOT SURE IF THIS IS BETTER THAN JUST USING
                %-2*epsilon_Z, in  the matrix, but this way we preserve
                %symmetry
                
                   
              %IN THIS VERSION WE ALWAYS HAVE DIRICHLET BC'S for top
              %electrode--> just include that here....
                %for Dirichet BC's, rhs should just = V_topBC --> we just
                %reset the bV(index) for these
%                 if  ( (i == floor((N+1)/2) || (i == floor((N+1)/2) + 1)) && (j == floor((N+1)/2) || (j == floor((N+1)/2) + 1)))
                    
                    bV(index) =  V_topBC(i,j); %note: this corresponds to BC eqn which just specifies the Dirichlet BC
                    %NOTE: use the V_topBC matrix to specify which i,j have a BC and which are 0 --> i.e. NOT really 0, but corresponds to Neumann BC's
                    %NOTE: V_topBC is non-zero only  for the i,j
                    %corresponding to the tip--> Dirichlet BC's
                    
                    %if not in tip region, are in insulating region, and no BC
                    %needs to be applied from rhs, b/c already taken care of in
                    %the matrix.
%                 end
                
            end
        end
    end
end
                        
                   
%Notice how much simpler rhs is with this ordering
%and including PBC's in the matrix!