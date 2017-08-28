function  bn = Setbn(B, n_initial, num_cell, num_elements, Cn, Un)

bn = zeros(num_elements,1);

%enforce boundary conditions through bp
bn(1) = 0;       %enforce, left size n is 0
bn(num_elements) = -B(1,num_elements+1)*n_initial;       

%introduce a net generation rate at 1 mesh point to the right of where
%introduce the bp generation rate
bn(ceil(num_cell/2.+1)) = -Cn*Un;