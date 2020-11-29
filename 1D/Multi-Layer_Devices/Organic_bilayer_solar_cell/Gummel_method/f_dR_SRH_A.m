function f_dR_SRH_A =  f_dR_SRH_A(E, n_full, p_full)
global gap_A q kb T T_tA cap_n cap_p N_LUMO N_HOMO Vt l N H_A

n1A = N_LUMO*exp(E/(kb*T));
p1A =	N_HOMO*exp((-gap_A*q-E)/(kb*T));
f_dR_SRH_A = H_A/(kb*T_tA)*exp(E/(kb*T_tA))*cap_n*cap_p*(n_full(l+1)*p_full(l)*N^2-N_LUMO*N_HOMO*exp(-gap_A/Vt))/(cap_n*(n_full(l+1)*N + n1A) + cap_p*(p_full(l)*N + p1A));
