%dBn(-dVi)/dVi
function dnegBni_dVi =  dnegBni_dVi(fullV)
global phi_c gap_elec Vt num_cell l

dV = zeros(1,num_cell+1);

for i = 2:num_cell+1                 %so dV(100) = dV(101: is at x=L) - dV(100, at x= l-dx).
    dV(i) = fullV(i)-fullV(i-1);     %these fullV's ARE ACTUALLY PSI PRIME; ALREADY = V/Vt, when solved poisson.
end

dV(num_cell+1) = dV(num_cell+1) + phi_c/Vt;          %injection step at right bndry (electrons electrode)
dV(l+1) = dV(l+1) + gap_elec/Vt;                     %energy step at interface

for i = 2:num_cell +1
    exp_term =  exp(dV(i));
    dnegBni_dVi(i) = exp_term*(-dV(i) + exp_term - 1)/(exp_term - 1)^2;
end