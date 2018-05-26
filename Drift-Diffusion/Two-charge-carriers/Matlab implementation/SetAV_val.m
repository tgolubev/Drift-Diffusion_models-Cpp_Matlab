function AV_val = SetAV_val(epsilon) 
global num_elements nx

AV = zeros(nx-2,3);

for i=1:num_elements     
    AV(i,2) = -2.*epsilon(i+1);     %i+1 b/c we fill matrix only corresponding to elements INSIDE the device
end

for i = 1:num_elements-1
    AV(i,1) = epsilon(i+1);    %1st element here corresponds to 2nd row of matrix
end
AV(num_elements, 1) = 0;     %last element is unused

for i = 2:num_elements
    AV(i,3) = epsilon(i+1);      %1st element here corresponds to the 1st row of matrix: so need to use i.e. epsilon corresponding to fullV(3) 
end


AV(1,3) = 0;                 %1st element of Ap(:,3) is unused

AV_val = spdiags(AV,-1:1,num_elements,num_elements);  %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.

