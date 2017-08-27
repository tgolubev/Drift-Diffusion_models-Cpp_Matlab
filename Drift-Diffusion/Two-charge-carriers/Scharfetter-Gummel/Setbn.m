function  bn = Setbn(B, n_initial, num_cell, num_elements, Cn, Un)

bn = zeros(num_elements,1);

%enforce boundary conditions through bp
bn(1) = -B(2,2)*n_initial;       %NOTE: scaled p_initial = 1 --> this is just like in fortran version
bn(num_elements) = 0;%-B(2,nx-1)*n(nx);       %ENFORCE RIGHT SIDE P IS = 0

%introduce a net generation rate at 1 mesh point to the right of where
%introduce the bp generation rate
bn(ceil(num_cell/2.+1)) = -Cn*Un;