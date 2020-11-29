function dFp_dV = dFp_dV(fullV, p_full, eps_V, Up)

global num_cell num_elements p_mob
%this will call Bernoulli, and recalculate Jn
%Will recalculate entire array--> using 1 eps_n value...

%compute fnc values to right and left of the node's V value...
B_p_plus = BernoulliFnc_p(fullV + eps_V/2.);
B_p_minus = BernoulliFnc_p(fullV - eps_V/2.);

% Calculate drift diffusion J's
% Use the SG definition--> but coeffs go into Cn and Cp...
for i = 1:num_cell        %verified that this is right
    Jp_plus(1,i+1) =  -p_mob*(p_full(i+1)*B_p_plus(1,i+1)-p_full(i)*B_p_plus(2,i+1));
    Jp_minus(1,i+1) =  -p_mob*(p_full(i+1)*B_p_minus(1,i+1)-p_full(i)*B_p_minus(2,i+1));
end

for i = 1:num_elements %so now this is recalculating continuity function, using the shifted V values
    dFp_dV(i) = (Continuity_p_fnc(Jp_plus(1,i+2), Jp_plus(1,i+1), Up(i)) - Continuity_p_fnc(Jp_minus(1,i+2), Jp_minus(1,i+1), Up(i)))/eps_V;
end
