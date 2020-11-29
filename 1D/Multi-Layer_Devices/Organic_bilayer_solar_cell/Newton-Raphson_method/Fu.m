function Fu = Net_Genfnc(U, k_ppd, zeta, Rtotal)

%this is treated as a function of U--> from here, the new U is found..

global l num_elements
%NOTE: b/c U is only non-zero at 2 points: l and l+1, when count real
%device indices, Fu has non-zero values only at 2 pts: Fu uses interior
%device points: so will be non-zero at l-1 and l (corresponding to l and
%l+1).

Fu = zeros(1,num_elements);

Fu(l-1) = U(l-1) - (k_ppd*zeta - Rtotal);  %note: k_ppd, zeta and Rtotal are not arrays
Fu(l) = U(l) - (k_ppd*zeta - Rtotal);