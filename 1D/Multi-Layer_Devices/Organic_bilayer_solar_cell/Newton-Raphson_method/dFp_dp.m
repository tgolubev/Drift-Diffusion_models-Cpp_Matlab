function dFp_dp = dFp_dp(B_p, p_full, Up, eps_p)

global num_cell p_mob num_elements
%here don't need to recalculate Bernoulli's b/c they are independent of n

for i = 1:num_cell        %verified that this is right
    %now plus and minus is determined by shifting n_full values...
    Jp_plus(1,i+1) =  p_mob*((p_full(i+1) + eps_p/2.)*B_p(1,i+1)-(p_full(i)+ eps_p/2.)*B_p(2,i+1));
    Jp_minus(1,i+1) = p_mob*((p_full(i+1) - eps_p/2.)*B_p(1,i+1)-(p_full(i)- eps_p/2.)*B_p(2,i+1));
end

for i = 1:num_elements %so now this is recalculating continuity function, using the shifted V values
    dFp_dp(i) = (Continuity_p_fnc(Jp_plus(1,i+2), Jp_plus(1,i+1), Up(i)) - Continuity_p_fnc(Jp_minus(1,i+2), Jp_minus(1,i+1), Up(i)))/eps_p;
end