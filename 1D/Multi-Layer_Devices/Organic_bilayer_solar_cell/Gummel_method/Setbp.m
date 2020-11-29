function bp = Setbp(B_p, p_full, Up)
global l Cp num_cell num_elements p_mob;

bp = zeros(num_elements,1);

%enforce boundary conditions through bp
bp(1,1) = -p_mob*B_p(1,2)*p_full(1);       
bp(num_elements, 1) = -p_mob*B_p(2,num_cell+1)*p_full(num_cell+1);       %ENFORCE RIGHT SIDE P IS = 0

%introduce a net generation rate somewhere in the middle
bp(l-1) = -Cp*Up;  