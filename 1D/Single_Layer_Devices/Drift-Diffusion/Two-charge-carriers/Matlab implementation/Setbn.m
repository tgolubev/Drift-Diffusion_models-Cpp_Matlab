function  bn = Setbn(B_n, n_full, Un)
global Cn num_cell num_elements n_mob;

bn = -Cn*Un;  

%enforce boundary conditions through bp
bn(1,1) = bn(1,1) - n_mob(2)*B_n(2,2)*n_full(1);       
bn(num_elements,1) = bn(num_elements,1) - n_mob(num_cell +1)*B_n(1,num_cell+1)*n_full(num_cell+1);      


