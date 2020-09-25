% RHS of matrix equation for hole charge density
function bp = Setbp(B_p, p_full, Up)
global Cp num_cell num_elements p_mob;

bp = -Cp*Up;  

%enforce boundary conditions for hole charge densities through bp
bp(1,1) = bp(1,1) - p_mob(2)*B_p(1,2)*p_full(1);       
bp(num_elements, 1) = bp(num_elements,1) - p_mob(num_cell+1)*B_p(2,num_cell+1)*p_full(num_cell+1);      

