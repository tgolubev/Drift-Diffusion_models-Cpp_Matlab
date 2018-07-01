% Setup of 2D electron continuity matrix using loop and spdiag
%Note: this is almost the same as Ap setup, except the Bn_pos and Bn_neg
% are switched everywhere

%NOTE MOBS INDICES WILL NEED TO BE CORRECTEED WHEN USE VARYING MOBILITIES

function An = SetAn_2D(n_mob, Bernoulli_n_values)

global num_elements N

%extract variables from struct, for brevity in eqns
Bn_posX = Bernoulli_n_values.Bn_posX;
Bn_negX = Bernoulli_n_values.Bn_negX;
Bn_posZ = Bernoulli_n_values.Bn_posZ;
Bn_negZ = Bernoulli_n_values.Bn_negZ;

An_val = zeros(num_elements, 5);   %this is a matrix which will just store the non-zero diagonals of 2D hole continuity eqn

%NOTE: index is not neccesarily equal to i (the x index of V), it is the
%index of the diagonal arrays.

%Lowest diagonal: corresponds to V(i, j-1)
i = 1;
j = 2;
for index = 1:N*(N-1)      %(1st element corresponds to Nth row  (number of elements = N*(N-1)
    An_val(index,1) = -((n_mob(i+1,j+1) + n_mob(i+1+1, j+1))/2.)*Bn_negZ(i+1,j+1);
    
     i = i+1;
    if (i > N)
        i = 1;  %reset i when reach end of subblock
        j = j+1;   %increment j
    end
end

%main lower diagonal (below main diagonal): corresponds to V(i-1,j)
i = 2;
j = 1;
for index = 1:num_elements-1      %(1st element corresponds to 2nd row)%NOTE: this is tricky!-->some elements are 0 (at the corners of the

%     if(mod(index,N) == 0)
%         An_val(index,2) = 0;   %these are the elements at subblock corners
%     else
    if (i > 1)
        An_val(index,2) = -((n_mob(i+1,j+1) + n_mob(i+1,j+1+1))/2.)*Bn_negX(i+1,j+1);
    end
    
    i = i+1;
    if (i > N)
        i = 1;
        j = j+1;
    end

end

%main diagonal: corresponds to V(i,j)
i = 1;
j = 1;
for index =  1:num_elements      
    An_val(index,3) = ((n_mob(i+1,j+1) + n_mob(i+1,j+1+1))/2.)*Bn_posX(i+1,j+1) + ((n_mob(i+1+1,j+1) + n_mob(i+1+1,j+1+1))/2.)*Bn_negX(i+1+1,j+1) + ((n_mob(i+1,j+1) + n_mob(i+1+1,j+1))/2.)*Bn_posZ(i+1,j+1) + ((n_mob(i+1,j+1+1) + n_mob(i+1+1,j+1+1))/2.)*Bn_negZ(i+1,j+1+1);
    
    i = i+1;
    if (i > N)
        i = 1;
        j = j+1;
    end
end

%main uppper diagonal: corresponds to V(i+1,j)
i = 1;
j = 1;
for index = 2:num_elements     %matlab fills this from the bottom (so i = 2 corresponds to 1st row in matrix)
%      
%     if(i == 0)     %i = 0, corresponds to this--> is actually  i = N
%         An_val(index,4) = 0;
%     else
    if (i > 0)
        An_val(index,4) = -((n_mob(i+1+1,j+1) + n_mob(i+1+1,j+1+1))/2.)*Bn_posX(i+1+1,j+1);
    end
    
     i=i+1;
    if (i > N-1)  %here are only N-1 elements per block and i starts from 1
        i = 0;
        j  = j+1;
    end
end

%far upper diagonal: corresponds to V(i,j+1)
i = 1;
j = 1;
for index = 1+N:num_elements      %matlab fills from bottom, so this starts at 1+N (since 1st element is in the 2nd subblock of matrix)
    An_val(index,5) = -((n_mob(i+1,j+1+1) + n_mob(i+1+1,j+1+1))/2.)*Bn_posZ(i+1,j+1+1);            %1st element corresponds to 1st row.   this has N^2 -N elements

    i = i+1;
    
    if (i > N)
        i = 1;
        j = j+1;
    end
end

%all not specified elements will remain zero, as they were initialized
%above.


An = spdiags(An_val, [-N -1 0 1 N], num_elements, num_elements); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
%diagonals  [-N -1 0 1 N].  -N and N b/c the far diagonals are in next
%subblocks, N diagonals away from main diag.
