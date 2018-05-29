% Setup of 2D Poisson matrix using loop and spdiag

function AV = SetAV_2D(epsilon)

global num_elements N

AV_val = zeros(num_elements, 5);   %this is a matrix which will just store the non-zero diagonals of 2D Poisson matrix

%These are listed in order from lowest diagonal to highest diagonal

%NOTE: index is not neccesarily equal to i (the x index of V), it is the
%index of the diagonal arrays.

%Lowest diagonal: corresponds to V(i, j-1)
for index = 1:N*(N-1)      % (1st element corresponds to Nth row  (number of elements = N*(N-1)
    i = mod(index,N);
    if(i ==0)                %the multiples of N correspond to last index
        i = N;
    end
    j = 1 + floor((index-1)/N);    %this is the y index of V which element corresponds to. 1+ floor((index-1)/N)determines which subblock this corresponds to and thus determines j, since the j's for each subblock are all the same.
    
    AV_val(index,1) = -(epsilon(i,j) + epsilon(i+1, j))/2.;
end

%NOTE: this is tricky!-->some elements are 0 (at the corners of the
%subblocks)
for index = 1:num_elements-1      %this is the lower diagonal (below main diagonal) (1st element corresponds to 2nd row)
    i = mod(index,N);         %this is x index of V which element corresponds to (note if this = 0, means these are the elements which are 0);
    j = 1 + floor((index-1)/N);
    
    if(mod(index, N) == 0)
        AV_val(index,2) = 0;   %these are the elements at subblock corners
    else
        AV_val(index,2) = -(epsilon(i,j) + epsilon(i,j+1))/2.;
    end
end

for index =  1:num_elements      %main diagonal
    i = mod(index,N);
    if(i ==0)                %the multiples of N correspond to last index
        i = N;
    end
    j = 1 + floor((index-1)/N);
    
    AV_val(index,3) = (epsilon(i+1,j) + epsilon(i+1,j+1))/2. + (epsilon(i,j) + epsilon(i,j+1))/2. + (epsilon(i,j+1) + epsilon(i+1,j+1))/2. + (epsilon(i,j) + epsilon(i+1,j))/2.;
end

for index = 2:num_elements      %main uppper diagonal, matlab fills this from the bottom (so i = 2 corresponds to 1st row in matrix)
    i = mod(index,N);
    j = 1 + floor((index-1)/N);
    
    if(mod(index-1,N) ==0)       %i-1 b/c indexed from 2
        AV_val(index,4) = 0;
    else
        AV_val(index,4) = -(epsilon(i+1,j) + epsilon(i+1,j+1))/2.;
    end
end

for index = 1+N:num_elements      %far upper diagonal, matlab fills from bottom, so this starts at 1+N (since 1st element is in the 2nd subblock of matrix)
    i = mod(index,N);
    if(i ==0)                %the multiples of N correspond to last index
        i = N;
    end
    j = 1 + floor((index-N)/N);
    
    AV_val(index,5) = -(epsilon(i,j+1) + epsilon(i+1,j+1))/2.;            %1st element corresponds to 1st row.   this has N^2 -N elements
end

%all not specified elements will remain zero, as they were initialized
%above.


AV = spdiags(AV_val, [-N -1 0 1 N], num_elements, num_elements); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
%diagonals  [-N -1 0 1 N].  -N and N b/c the far diagonals are in next
%subblocks, N diagonals away from main diag.

