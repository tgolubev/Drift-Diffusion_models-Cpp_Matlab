function R_Langevin_array = R_Langevin(n_full, p_full)
global k_rec n1 p1 num_cell N

R_Langevin_array = zeros(num_cell,1);

for i = 2:num_cell      
    R_Langevin_array(i,1) = k_rec*(N*n_full(i)*N*p_full(i) - n1*p1);
    
    %negative recombo rates are not physical
    if R_Langevin_array(i) < 0 
        R_Langevin_array(i) = 0; 
    end
end

R_Langevin_array = R_Langevin_array(2:num_cell);  %to make sure aligned with G
    