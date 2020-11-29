function FnT_element = nTraps_fnc(ni, nTi)

global H_A T T_tA N N_LUMO 

FnT_element = nTi - (H_A/N)*(ni*N/N_LUMO).^(T/T_tA);  %will be array of num_elements length

