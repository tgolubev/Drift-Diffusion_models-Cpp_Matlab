function R = R(n_full, p_full)
global gap_A gap_D q kb T N N_LUMO N_HOMO l;

%Calcurate R_SRH with acceptor trap
a = -gap_A*q;
b = kb*T*log(n_full(l+1)*N/N_LUMO);  %the log here is an ln
n = 100;       %I verified that 100 is enough to give reasonable result for the integration 
sumk = 0.0;

for k= 1:n-1
sumk = sumk + f_dR_SRH_A(a + k*(b-a)/n, n_full, p_full);    %the first term here is E in f_dR_SRH_A
end
R_SRH_A = (b-a)/n*(f_dR_SRH_A(a, n_full, p_full)/2.0D0 + sumk + f_dR_SRH_A(b, n_full, p_full)/2.0D0);

%Calcurate R_SRH with donor trap
a = kb*T*log(p_full(l)*N/N_HOMO);
b = -gap_D*q;
n = 100;
sumk = 0.0D0;

for k=1:n-1
sumk = sumk + f_dR_SRH_D(a + k*(b-a)/n, n_full, p_full);
end
R_SRH_D = (b-a)/n*(f_dR_SRH_D(a, n_full, p_full)/2.0D0 + sumk + f_dR_SRH_D(b, n_full, p_full)/2.0D0);

R_SRH = R_SRH_A + R_SRH_D;
R = R_SRH;