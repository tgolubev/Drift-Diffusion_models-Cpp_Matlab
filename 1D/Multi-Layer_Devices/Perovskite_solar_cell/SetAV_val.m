function AV_val = SetAV_val(epsilon) 
global num_elements nx l_HTL_int l_ETL_int l_BCP_int

AV = zeros(nx-2,3);

for i=1:num_elements     
    AV(i,2) = -2.*epsilon(i+1);    % i+1 b/c we fill matrix only corresponding to elements INSIDE the device
end

for i = 1:num_elements-1
    AV(i,1) = epsilon(i+1);  % 1st element here corresponds to 2nd row of matrix
end
AV(num_elements, 1) = 0;     % last element is unused


for i = 2:num_elements
    AV(i,3) = epsilon(i+1);      %1st element here corresponds to the 1st row of matrix
end

% Rows in the matrix corresponding to interfaces between layers need
% special treatement. When we discretize over the interface, 
% we need to use the appropriate epsilon value
 AV(l_HTL_int,3) = (epsilon(1)); 
 AV(l_HTL_int,2) = -(epsilon(floor(num_elements/2)) + epsilon(1));  

 AV(l_ETL_int,3) = (epsilon(floor(num_elements/2)));
 AV(l_ETL_int ,2) = -(epsilon(l_ETL_int+1) + epsilon(floor(num_elements/2)));
 
 AV(l_BCP_int,3) = (epsilon(l_ETL_int+1)); 
 AV(l_BCP_int,2) = -(epsilon(num_elements) + epsilon(l_ETL_int +1));

AV(1,3) = 0;                 % 1st element of Ap(:,3) is unused

AV_val = spdiags(AV,-1:1,num_elements,num_elements);  