function Bp = BernoulliFnc_p(fullV)

global num_cell 

dV = zeros(num_cell+1);
Bp = zeros(2, num_cell+1);

%Bernoulli fnc.

 for i = 2:num_cell+1                 
     dV(i) = fullV(i)-fullV(i-1);     
 end
 
 for i = 2:num_cell+1
     Bp(1,i) = dV(i)/(exp(dV(i))-1.0);    %B(+dV)
     Bp(2,i) = Bp(1,i)*exp(dV(i));          %B(-dV)
 end
