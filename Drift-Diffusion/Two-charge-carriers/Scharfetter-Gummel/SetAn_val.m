function An_val = SetAn_val(num_elements, B)

An =  zeros(num_elements,3);
                
%BE CAREFUL!: the setup of the  sparse matrix elements is
%non-trivial: last entry  of Ap(:,1) 1st entry of
%Ap(:,3) are unused b/c off diagonals have 1 less element than main
%diagonal!

for i=1:num_elements     %I verified that ordering of columns is correct!
    An(i,2) = -(B(1,i+1) + B(2,i+1+1));     %I HAVE VERIFIED THAT ALL THE dV's properly  match up with  indices! All the extra +1's are b/c B starts at i=2 (1st pt inside device).
end
for i = 1:num_elements-1
    An(i,1) = B(2,i+1+1);    %should be i+1+1 b/c have B(dVi)*p_i-1  (see fortran reindexed version)
end
An(num_elements, 1) = 0;     %last element is unused
for i = 2:num_elements
    An(i,3) = B(1,i+1);      %1st element here corresponds to the 1st row of matrix: so need to use B3 --> corresponding to p(3)
end
An(1,3) = 0;                 %1st element of Ap(:,3) is unused


An_val = spdiags(An,-1:1,num_elements,num_elements);