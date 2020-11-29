% This is the Bernoulli function for electrons from Scharfetter-Gummel discretization
function B_n = BernoulliFnc_n(fullV)
global num_cell Vt ETL_int_CBstep l_ETL_int HTL_int_CBstep l_HTL_int BCP_int_CBstep l_BCP_int

dV = zeros(1,num_cell+1);
B_n = zeros(2, num_cell+1);

 for i = 2:num_cell+1                 
     dV(i) = fullV(i)-fullV(i-1);   
 end

% energy barriers
% injection barriers at contacts should not be incuded b/c are already in
% the BCs
dV(l_HTL_int +1) = dV(l_HTL_int +1)  + HTL_int_CBstep/Vt;    
dV(l_BCP_int+1) = dV(l_BCP_int+1) + BCP_int_CBstep/Vt;  
dV(l_ETL_int +1) = dV(l_ETL_int +1)  +  ETL_int_CBstep/Vt;  

 for i = 2:num_cell+1
     B_n(1,i) = dV(i)/(exp(dV(i))-1.0);    % B(+dV)
     B_n(2,i) = B_n(1,i)*exp(dV(i));       % B(-dV)
 end
