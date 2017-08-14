function B = BernoulliFnc(nx, fullV, Vt)

dV = zeros(nx-1);
B = zeros(2, nx-1);

%Bernoulli fnc.

%NOTE: THESE DEFINITIONS ARE DIFFERENT THAN IN DD-BI code


 for i = 1:nx-1               %so dV(100) = dV(101: is at x=L) - dV(100, at x= l-dx).
    dV(i) = fullV(i+1)-fullV(i);
    B(1,i) = dV(i)/(exp(dV(i))-1.0);    %B(+dV)
    B(2,i) = B(1)*exp(dV(i));          %B(-dV)
 end
