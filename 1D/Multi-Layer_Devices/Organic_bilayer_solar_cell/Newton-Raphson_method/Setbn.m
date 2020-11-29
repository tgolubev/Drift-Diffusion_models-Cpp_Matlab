function  bn = Setbn(B_n, n_full, Un)
global l Cn_gummel num_cell num_elements n_mob;

bn = zeros(num_elements,1);

%enforce boundary conditions through bp
bn(1) = -n_mob*B_n(2,2)*n_full(1);       
bn(num_elements) = -n_mob*B_n(1,num_cell+1)*n_full(num_cell+1);       %THERE WAS A MISTAKE HERE!! HAD num_elements +1 instead of num_cell+1! so weren't using the electron injection barrier! b/c num_elements +1 = num_cell

%introduce a net generation rate at 1 mesh point to the right of where
%introduce the bp generation rate
bn(l) = -Cn_gummel*Un;  
