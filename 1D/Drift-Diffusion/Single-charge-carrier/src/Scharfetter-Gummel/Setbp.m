function bp = Setbp(B, p_initial, num_cell, num_elements, Cp, U)

bp = zeros(num_elements,1);

%enforce boundary conditions through bp
bp(1) = -B(1,2)*p_initial;       %NOTE: scaled p_initial = 1 --> this is just like in fortran version
bp(num_elements) = 0;%-B(1,nx-1)*p(nx);       %ENFORCE RIGHT SIDE P IS = 0

%introduce a net generation rate somewhere in the middle
bp(ceil(num_cell/2.)) = -Cp*U;