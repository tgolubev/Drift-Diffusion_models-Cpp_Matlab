function R_Langevin_array = R_Langevin(n_full, p_full)
global k_rec n1 p1 l_HTL_int l_ETL_int num_cell N

R_Langevin_array = zeros(num_cell+1,1);

for i = l_HTL_int +1:l_ETL_int       %for within pervoskite layer: layer starts at l_HTL_+1....
    R_Langevin_array(i) = k_rec*(N*n_full(i)*N*p_full(i) - n1*p1);
end

R_Langevin_array = R_Langevin_array(2:num_cell);  % to make sure aligned with G
    
