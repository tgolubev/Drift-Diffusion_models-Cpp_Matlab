%dBp(-dVi)/dVi
function dnegBpi_dVi =  dnegBpi_dVi(fullV)
global phi_a gap_hole Vt num_cell l

dV = zeros(1,num_cell+1);

for i = 2:num_cell+1                 %so dV(100) = dV(101: is at x=L) - dV(100, at x= l-dx).
    dV(i) = fullV(i)-fullV(i-1);     %these fullV's = V/Vt, when solved poisson.
end

dV(2) = dV(2) + phi_a/Vt;               %injection step at left bndry (holes electrode)
dV(l+1) = dV(l+1) + gap_hole/Vt;

for i = 2:num_cell +1
    exp_term =  exp(dV(i));
    dnegBpi_dVi(i) = exp_term*(-dV(i) + exp_term - 1)/(exp_term - 1)^2;
end