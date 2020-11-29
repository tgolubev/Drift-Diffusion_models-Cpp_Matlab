function B_p = BernoulliFnc_p(fullV)
global num_cell phi_a Vt l gap_hole

dV = zeros(1,num_cell+1);
B_p = zeros(2, num_cell+1);

%Bernoulli fnc.

 for i = 2:num_cell+1                 %so dV(100) = dV(101: is at x=L) - dV(100, at x= l-dx).
     dV(i) = fullV(i)-fullV(i-1);     %these fullV's ARE ACTUALLY PSI PRIME; ALREADY = V/Vt, when solved poisson.
 end
 
 dV(2) = dV(2) + phi_a/Vt;               %injection step at left bndry (holes electrode) 
 %dV(l+1) = dV(l+1) + gap_hole/Vt;        %energy step at interface

 dV(l+1) = dV(l+1) + gap_hole/Vt;  
 
 for i = 2:num_cell+1
     B_p(1,i) = dV(i)/(exp(dV(i))-1.0);    %B(+dV)
     %B(2,i) = -dV(i)/(exp(-dV(i))-1.0); %THIS IS EQUIVALENT  TO OTHE
     %BELOW EXPERSSION
     B_p(2,i) = B_p(1,i)*exp(dV(i));          %B(-dV)
 end
