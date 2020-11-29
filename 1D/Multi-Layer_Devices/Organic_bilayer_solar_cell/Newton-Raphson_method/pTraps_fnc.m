function FpT_element = pTraps_fnc(pi, pTi)

global H_D T T_tD N N_HOMO 

FpT_element = pTi - (H_D/N)*(pi*N/N_HOMO).^(T/T_tD);  %will be array of num_elements length