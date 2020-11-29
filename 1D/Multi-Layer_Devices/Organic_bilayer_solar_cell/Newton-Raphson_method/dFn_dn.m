function dFn_dn = dFn_dn(B_n, n_full, Un, eps_n)

global num_cell n_mob num_elements
%here don't need to recalculate Bernoulli's b/c they are independent of n

for i = 1:num_cell        %verified that this is right 
    %now plus and minus is determined by shifting n_full values...
    Jn_plus(1,i+1) =  n_mob*((n_full(i+1) + eps_n/2.)*B_n(1,i+1)-(n_full(i)+ eps_n/2.)*B_n(2,i+1)); 
    Jn_minus(1,i+1) = n_mob*((n_full(i+1) - eps_n/2.)*B_n(1,i+1)-(n_full(i)- eps_n/2.)*B_n(2,i+1)); 
end

for i = 1:num_elements %so now this is recalculating continuity function, using the shifted V values
    dFn_dn(i) = (Continuity_n_fnc(Jn_plus(1,i+2), Jn_plus(1,i+1), Un(i)) - Continuity_n_fnc(Jn_minus(1,i+2), Jn_minus(1,i+1), Un(i)))/eps_n;
end
