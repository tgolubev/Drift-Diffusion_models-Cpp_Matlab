function An_val = SetAn_val(B_n)
global num_elements n_mob 

An =  zeros(num_elements,3);
                
%last entry  of An(:,1) 1st entry of
%An(:,3) are unused b/c off diagonals have 1 less element than main
%diagonal!

for i=1:num_elements    
    An(i,2) = -(n_mob(i+1)*B_n(1,i+1) + n_mob(i+1+1)*B_n(2,i+1+1));    
end
for i = 1:num_elements-1
    An(i,1) = n_mob(i+1+1)*B_n(2,i+1+1);    %should be i+1+1 b/c have B(dVi)*p_i-1  
end
An(num_elements, 1) = 0;     %last element is unused
for i = 2:num_elements
    An(i,3) = n_mob(i+1)*B_n(1,i+1);      %1st element here corresponds to the 1st row of matrix: so need to use B3 --> corresponding to p(3)
end
An(1,3) = 0;                 %1st element of An(:,3) is unused


An_val = spdiags(An,-1:1,num_elements,num_elements);