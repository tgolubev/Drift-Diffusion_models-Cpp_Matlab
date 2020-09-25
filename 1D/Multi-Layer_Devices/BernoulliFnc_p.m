% This is the Bernoulli function for holes from Scharfetter-Gummel discretization
function B_p = BernoulliFnc_p(fullV)
global num_cell l_HTL_int Vt HTL_int_VBstep phi_a l_ETL_int ETL_int_VBstep phi_c L_int L_int_eachside BCP_int_VBstep l_BCP_int

dV = zeros(1,num_cell+1);
B_p = zeros(2, num_cell+1);

 for i = 2:num_cell+1                 
     dV(i) = fullV(i)-fullV(i-1);  
 end

% energy barriers
% injection barriers at contacts should not be incuded b/c are already in
% the BCs
dV(l_HTL_int +1 ) = dV(l_HTL_int +1) + HTL_int_VBstep/Vt; 
dV(l_ETL_int +1) = dV(l_ETL_int +1)  + ETL_int_VBstep/Vt; 
dV(l_BCP_int+1) = dV(l_BCP_int+1) + BCP_int_VBstep/Vt;
 
 for i = 2:num_cell+1
     B_p(1,i) = dV(i)/(exp(dV(i))-1.0);    % B(+dV)
     B_p(2,i) = B_p(1,i)*exp(dV(i));       % B(-dV)
 end
