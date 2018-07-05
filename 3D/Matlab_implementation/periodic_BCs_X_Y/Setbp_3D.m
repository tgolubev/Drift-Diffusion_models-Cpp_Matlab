function bp = Setbp_3D(Bernoulli_p_values, p_mob, Up)
global Cp num_elements p_topBC p_bottomBC N;

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
for i = 1:N+1
        for j = 1:N+1
            for k = 1:N
                index = index + 1;
                if (k == 1)  %bottom BC                  
                    bp(index) = bp(index) + p_mob(i+1,j+1,0+1)*p_bottomBC(i,j)*Bp_posZ(i+1,j+1,k+1);
                end
                %middle elements have no BC's
                if (k == N) %top BC
                    bp(index) = bp(index) + p_mob(i+1,j+1,N+1)*p_topBC(i,j)*Bp_negZ(i+1,j+1,k+1+1);
                end
            end
        end
end