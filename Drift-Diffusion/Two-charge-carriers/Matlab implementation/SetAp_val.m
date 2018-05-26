function Ap_val = SetAp_val(B_p) 
global num_elements p_mob 

Ap =  zeros(num_elements,3);
                
% last entry  of Ap(:,1) 1st entry of
%Ap(:,3) are unused b/c off diagonals have 1 less element than main
%diagonal!

%NOTE: Mobilities indices just correspond to the the dV indices which are
%same as B indices

for i=1:num_elements     
    Ap(i,2) = -(p_mob(i+1)*B_p(2,i+1) + p_mob(i+1+1)*B_p(1,i+1+1));    
end

for i = 1:num_elements-1
    Ap(i,1) = p_mob(i+1+1)*B_p(1,i+1+1);    %should be i+1+1 b/c have B(dVi)*p_i-1  
end
Ap(num_elements, 1) = 0;     %last element is unused
for i = 2:num_elements
    Ap(i,3) = p_mob(i+1)*B_p(2,i+1);      %1st element here corresponds to the 1st row of matrix: so need to use B3 --> corresponding to p(3)--> I mean fullp(3) here
end

Ap(1,3) = 0;                 %1st element of Ap(:,3) is unused


Ap_val = spdiags(Ap,-1:1,num_elements,num_elements);
