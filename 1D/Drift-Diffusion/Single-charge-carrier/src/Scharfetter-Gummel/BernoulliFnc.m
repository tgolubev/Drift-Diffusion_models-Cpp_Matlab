function B = BernoulliFnc(num_cell, fullV)

dV = zeros(1,num_cell+1);
B = zeros(2, num_cell+1);

%Bernoulli fnc.

 for i = 2:num_cell+1                 %so dV(100) = dV(101: is at x=L) - dV(100, at x= l-dx).
     dV(i) = fullV(i)-fullV(i-1);     %these fullV's ARE ACTUALLY PSI PRIME; ALREADY = V/Vt, when solved poisson.
 end
 
 for i = 2:num_cell+1
     B(1,i) = dV(i)/(exp(dV(i))-1.0);    %B(+dV)
     %B(2,i) = -dV(i)/(exp(-dV(i))-1.0); %THIS IS EQUIVALENT  TO OTHE
     %BELOW EXPERSSION
     B(2,i) = B(1,i)*exp(dV(i));          %B(-dV)
 end
