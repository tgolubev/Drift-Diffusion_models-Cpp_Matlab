function Fv_element = Poissonfnc(fullV, ni, pi,nTi,  pTi, i)
%NOTE: define imput argument n(i) = ni --> is more robust...
%since only need to pass 1 n and p value, don't pass the  whole array!
%NOTE: need to pass it i, so it knows what element are computing!

%NOTE: need fullV, b/c need boundary fullV values when do laplacian

global epsilon q num_elements CV
%NOTE: epsilon = (relative permitivity)*epsilon_0


%the function will be defined for INSIDE the device: from cell 2 to
%num_cell --> or since V, n, and, p are already for INSIDE device, for all
%those elements...

% eps_r = epsilon/epsilon_0;


Fv_element = ((fullV(i) -2*fullV(i+1) + fullV(i+2)) - CV*((ni+nTi)-(pi+pTi)));         %the fullV values use +1 everywhere, b/c starts from 1 = endpoint, whereas for n and p, 1 = 1st interior point
