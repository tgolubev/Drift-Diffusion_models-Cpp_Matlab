function dFn_dV = dFn_dV(fullV, n_full, eps_V, Un)

global num_cell num_elements n_mob
%this will call Bernoulli, and recalculate Jn
%Will recalculate entire array--> using 1 eps_n value...

%compute fnc values to right and left of the node's V value...
B_n_plus = BernoulliFnc_n(fullV + eps_V/2.);
B_n_minus = BernoulliFnc_n(fullV - eps_V/2.);

% Calculate drift diffusion J's
% Use the SG definition--> BUT COEFFS go into Cn and Cp...
for i = 1:num_cell        %verified that this is right
    Jn_plus(1,i+1) =  n_mob*(n_full(i+1)*B_n_plus(1,i+1)-n_full(i)*B_n_plus(2,i+1));
    Jn_minus(1,i+1) =  n_mob*(n_full(i+1)*B_n_minus(1,i+1)-n_full(i)*B_n_minus(2,i+1));
end

for i = 1:num_elements %so now this is recalculating continuity function, using the shifted V values
    dFn_dV(i) = (Continuity_n_fnc(Jn_plus(1,i+2), Jn_plus(1,i+1), Un(i)) - Continuity_n_fnc(Jn_minus(1,i+2), Jn_minus(1,i+1), Un(i)))/eps_V;
end




