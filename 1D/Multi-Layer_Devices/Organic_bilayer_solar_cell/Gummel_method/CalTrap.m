function [pT_full, nT_full] = CalTrap(p_full, n_full)  %allows to return 2 values with this notation
global H_D H_A N N_HOMO N_LUMO T T_tD T_tA;

pT_full = H_D*(p_full*N/N_HOMO).^(T/T_tD);   %USE .^ for element wise (scalar) power
nT_full = H_A*(n_full*N/N_LUMO).^(T/T_tA);

pT_full = pT_full/N;
nT_full = nT_full/N;

