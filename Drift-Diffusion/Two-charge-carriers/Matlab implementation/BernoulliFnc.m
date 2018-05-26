function B = BernoulliFnc(fullV)

global num_cell

dV = zeros(num_cell+1);
B = zeros(2, num_cell+1);

%Bernoulli fnc.

 for i = 2:num_cell+1                 
     dV(i) = fullV(i)-fullV(i-1);     
 end
 
 for i = 2:num_cell+1
     B(1,i) = dV(i)/(exp(dV(i))-1.0);    %B(+dV)
     B(2,i) = B(1,i)*exp(dV(i));          %B(-dV)
 end
