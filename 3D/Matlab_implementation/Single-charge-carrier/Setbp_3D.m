function bp = Setbp_3D(Bernoulli_p_values, p_mob, Up)
global Cp num_elements p_topBC p_bottomBC Nx Ny Nz;

% bp = zeros(num_elements,1);
%Note: assuming for BC's that mobility  on bndry is same as mobility just
%inside the boundary (so not doing any mobs averaging there)

%extract variables from struct (for brevity in eqns)
% Bp_posX = Bernoulli_p_values.Bp_posX;
% Bp_negX = Bernoulli_p_values.Bp_negX;
% Bp_posY = Bernoulli_p_values.Bp_posY;
% Bp_negY = Bernoulli_p_values.Bp_negY;
Bp_posZ = Bernoulli_p_values.Bp_posZ;
Bp_negZ = Bernoulli_p_values.Bp_negZ;

%calculate main part here
bp = Cp*Up;

%add on BC's
index = 0;
for i = 1:Nx+1
    for j = 1:Ny+1
        for k = 1:Nz+1
            index = index + 1;
            if (k == 1)  %bottom BC
                bp(index) = bp(index) + p_mob(i+1,j+1,0+1)*p_bottomBC(i,j)*Bp_posZ(i+1,j+1,k+1);
            end
            %middle elements have no BC's
            if (k == Nz+1) %top BC
                
%                 bp(index) = bp(index)/2;   %this is due to I used -p_mob_Z_avg.*Bp_posZ instead of -2*p_mob_Z_avg.*Bp_posZ in Poisson for neuman BC
                %NOT SURE IF THIS IS BETTER THAN JUST USING
                %-2*p_mob_Z_avg.*Bp_posZ, in  the matrix, ESPECIALLY
                %SINCE Ap matrix is not symmetry already--> doing this
                %doessn't preserve symmetry....
                
                
              %IN THIS VERSION WE ALWAYS HAVE DIRICHLET BC'S for top
              %electrode--> just include that here....
%                 if  ( (i == floor((N+1)/2) || (i == floor((N+1)/2) + 1)) && (j == floor((N+1)/2) || (j == floor((N+1)/2) + 1)))
                    bp(index) = p_topBC(i,j);  %note: this corresponds to BC eqn which just specifies the Dirichlet BC
                    
                    %NOTE: p_topBC is non-zero only  for the i,j
                    %corresponding to the tip--> Dirichlet BC's
                    
                    %if not in tip region, are in insulating region, and no BC
                    %needs to be applied from rhs, b/c already taken care of in
                    %the matrix.
%                 end
            end
            
            
        end
    end
end
end