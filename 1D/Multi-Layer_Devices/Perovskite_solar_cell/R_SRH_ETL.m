function R_SRH_ETL_array = R_SRH_ETL(n_full, p_full)
global ETL_traps num_elements l_ETL_int L_int_eachside cap_n cap_p  n1 p1 N L_int;

R_SRH_ETL_array = zeros(num_elements,1);

for i = l_ETL_int-9:l_ETL_int 
    R_SRH_ETL_array(i) = cap_n*cap_p*ETL_traps*(N*n_full(i)*N*p_full(i) - n1*p1)/(cap_n*(N*n_full(i) + n1) + (cap_p*(N*p_full(i) + p1)));
    if(R_SRH_ETL_array(i,1) < 0.0)
        R_SRH_ETL_array(i,1) = 0.0;  % negative recombination rate is not physical
    end
end

