function f_dR_SRH_D = f_dR_SRH_D(E, n_full, p_full)
global gap_D q kb T T_tD cap_n cap_p N_LUMO N_HOMO Vt l N H_D


n1D = N_LUMO*exp((-gap_D*q-E)/(kb*T));
p1D = N_HOMO*exp(E/(kb*T));
f_dR_SRH_D = -H_D/(kb*T_tD)*exp(E/(kb*T_tD))*cap_n*cap_p*(n_full(l+1)*p_full(l)*N^2 - N_LUMO*N_HOMO*exp(-gap_D/Vt))/(cap_n*(n_full(l+1)*N + n1D) + cap_p*(p_full(l)*N + p1D));