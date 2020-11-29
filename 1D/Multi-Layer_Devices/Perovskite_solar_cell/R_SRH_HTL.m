function R_SRH_HTL_array = R_SRH_HTL(n_full, p_full)
global HTL_traps num_elements l_HTL_int  cap_n cap_p n1 p1 N L_int L_int_eachside;

R_SRH_HTL_array = zeros(num_elements,1);

for i = l_HTL_int+1:l_HTL_int+10 
    R_SRH_HTL_array(i,1) = cap_n*cap_p*HTL_traps*(N*n_full(i).*N*p_full(i) - n1*p1)/(cap_n*(N*n_full(i) + n1) + (cap_p*(N*p_full(i) + p1)));
    if(R_SRH_HTL_array(i,1) < 0.0)
        R_SRH_HTL_array(i,1) = 0.0;  % negative recombination rate is not physical
    end
end


