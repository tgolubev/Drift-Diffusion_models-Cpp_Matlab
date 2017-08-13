function Ap_val = SetAp_val(num_cell, fullV, p, Vt) 

Ap =  zeros(num_cell,3);

%Bernoulli fnc.

 for i = 2:num_cell+1                 %num_pts including boundaries = num_cell+1
    dV(i) = fullV(i)-fullV(i-1);
    B(1,i) = dV(i)/(exp(dV(i))-1.0);    %B(+dV)
    B(2,i) = B(1)*exp(dV(i));          %B(-dV)
 end
                
 for i=2:num_cell
    Ap(i,1) = B(1,i);     %lower diagonal
    Ap(i,2) = B(2,i) + B(1,i+1);  %main diagonal
    Ap(i,3) = B(2,i+1);     %upper diagonal   CHECK MAKE SURE THIS IS RIGHT!
 end
 
Ap_val = spdiags(Ap,-1:1,num_cell-1,num_cell-1); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.