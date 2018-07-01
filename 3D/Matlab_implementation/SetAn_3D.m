% Setup of 3D electron continuity matrix using loop and spdiag
%Note: this is almost the same as Ap setup, except the Bn_pos and Bn_neg
% are switched everywhere

function An = SetAn_3D(n_mob, Bernoulli_n_values)

global num_elements N

%extract variables from struct, for brevity in eqns
Bn_posX = Bernoulli_n_values.Bn_posX;
Bn_negX = Bernoulli_n_values.Bn_negX;
Bn_posY = Bernoulli_n_values.Bn_posY;
Bn_negY = Bernoulli_n_values.Bn_negY;
Bn_posZ = Bernoulli_n_values.Bn_posZ;
Bn_negZ = Bernoulli_n_values.Bn_negZ;

An_val = zeros(num_elements, 7);   %this is a matrix which will just store the non-zero diagonals of 2D hole continuity eqn

%NOTE: index is not neccesarily equal to i (the x index of V), it is the
%index of the diagonal arrays.
%--------------------------------------------------------------------------
%Lowest diagonal: 
i = 1;
j = 1;
k = 2;  %these are the 1st indices OF MAIN DIAG (or rhs), that 1st element of lowest diag corresponds to.
for index = 1:N^3 - N^2      % (1st element corresponds to Nth row  (number of elements = N^3 - N^2) 
    An_val(index,1) = -((n_mob(i+1,j+1,k+1) + n_mob(i+1+1,j+1,k+1) + n_mob(i+1, j+1+1, k+1) + n_mob(i+1+1,j+1+1,k+1))/4.)*Bn_negZ(i+1,j+1,k+1);

    i = i+1;
    if (i > N)
        i = 1;  %reset i when reach end of subblock
        j = j+1;   
    end
     if (j > N)  
        j = 1;
        k = k+1;
    end
    
end  
%--------------------------------------------------------------------------
%lower diagonal
i = 1;
j = 2;  %NOTE: this are the rhs indices (= main diag indices) which these elements correspond to.
k = 1;
for index = 1:N^3 - N      
    if (j > 1)
        An_val(index,2) = -((n_mob(i+1,j+1, k+1) + n_mob(i+1+1, j+1, k+1) + n_mob(i+1, j+1, k+1+1) + n_mob(i+1+1, j+1, k+1+1))/4.)*Bn_negY(i+1,j+1,k+1);
    end
    
    i = i+1;
    if (i > N)
        i = 1;  %reset i when reach end of subblock
        j = j+1;
    end
    if (j > N)  
        j = 1;
        k = k+1;
    end
end      

%--------------------------------------------------------------------------
%main lower diagonal (below main diagonal)
i = 2;
j = 1;
k = 1;
for index = 1:N^3 - 1    
    if (i > 1)
        An_val(index,3) = -((n_mob(i+1,j+1, k+1) + n_mob(i+1,j+1+1, k+1) + n_mob(i+1,j+1,k+1+1) + n_mob(i+1,j+1+1,k+1+1))/4.)*Bn_negX(i+1,j+1,k+1);
    end
    
    i = i+1;
    if (i > N)
        i = 1;
        j = j+1;
    end
    if (j > N)
        j = 1;
        k = k+1;
    end
end
%--------------------------------------------------------------------------
%main diagonal
i = 1;
j = 1;
k = 1;
for index =  1:num_elements          
    An_val(index,4) = ((n_mob(i+1,j+1,k+1) + n_mob(i+1+1,j+1,k+1) + n_mob(i+1, j+1+1, k+1) + n_mob(i+1+1,j+1+1,k+1))/4.)*Bn_posZ(i+1,j+1,k+1) ...
                    + ((n_mob(i+1,j+1, k+1) + n_mob(i+1+1, j+1, k+1) + n_mob(i+1, j+1, k+1+1) + n_mob(i+1+1, j+1, k+1+1))/4.)*Bn_posY(i+1,j+1,k+1) ...
                    + ((n_mob(i+1,j+1, k+1) + n_mob(i+1,j+1+1, k+1) + n_mob(i+1,j+1,k+1+1) + n_mob(i+1,j+1+1,k+1+1))/4.)*Bn_posX(i+1,j+1,k+1) ...
                    + ((n_mob(i+1+1,j+1,k+1) + n_mob(i+1+1,j+1+1,k+1) + n_mob(i+1+1,j+1,k+1+1) + n_mob(i+1+1,j+1+1,k+1+1))/4.)*Bn_negX(i+1+1,j+1,k+1) ...
                    + ((n_mob(i+1, j+1+1, k+1) + n_mob(i+1+1, j+1+1, k+1) + n_mob(i+1, j+1+1, k+1+1) + n_mob(i+1+1, j+1+1, k+1+1))/4.)*Bn_negY(i+1,j+1+1,k+1) ...
                    + ((n_mob(i+1,j+1, k+1+1) + n_mob(i+1+1,j+1,k+1+1) + n_mob(i+1,j+1+1,k+1+1) + n_mob(i+1+1,j+1+1,k+1+1))/4.)*Bn_negZ(i+1,j+1,k+1+1);

    i = i+1;
    if (i > N)
        i = 1;
        j = j+1;
    end
    if (j > N)
        j = 1;
        k = k+1;
    end

end
%--------------------------------------------------------------------------
%main uppper diagonal
i = 1;
j = 1;
k = 1;
for index = 2:N^3 - 1 +1    %matlab fills this from the bottom (so i = 2 corresponds to 1st row in matrix)
    
    if (i > 0)
        An_val(index,5) = -((n_mob(i+1+1,j+1,k+1) + n_mob(i+1+1,j+1+1,k+1) + n_mob(i+1+1,j+1,k+1+1) + n_mob(i+1+1,j+1+1,k+1+1))/4.)*Bn_posX(i+1+1,j+1,k+1);
    end
    
    i=i+1;
    if (i > N-1)  %here are only N-1 elements per block and i starts from 1
        i = 0;
        j  = j+1;
    end
    if (j > N)
        j = 1;
        k = k+1;
    end
    
end      
%--------------------------------------------------------------------------
%upper diagonal
i = 1;
j = 1;
k = 1;
for index = 1+N:N^3 - N +N 
    if (j > 0)
        An_val(index, 6) = -((n_mob(i+1, j+1+1, k+1) + n_mob(i+1+1, j+1+1, k+1) + n_mob(i+1, j+1+1, k+1+1) + n_mob(i+1+1, j+1+1, k+1+1))/4.)*Bn_posY(i+1,j+1+1,k+1);
    end

     i = i+1;
     if (i > N)
        i = 1;
        j = j+1;
     end
     if (j > N-1)
        j = 0;
        k = k+1;
     end
end
%--------------------------------------------------------------------------
%far upper diagonal
i = 1;
j = 1;
k = 1;
for index = 1+N*N:num_elements      %matlab fills from bottom, so this starts at 1+N (since 1st element is in the 2nd subblock of matrix) 
    An_val(index,7) = -((n_mob(i+1,j+1, k+1+1) + n_mob(i+1+1,j+1,k+1+1) + n_mob(i+1,j+1+1,k+1+1) + n_mob(i+1+1,j+1+1,k+1+1))/4.)*Bn_posZ(i+1,j+1,k+1+1);         
    
    i = i+1;  
    if (i > N)
        i = 1;
        j = j+1;
    end
     if (j > N)
        j = 1;
        k = k+1;
    end
end 

%all not specified elements will remain zero, as they were initialized
%above.


An = spdiags(An_val, [-N^2 -N -1 0 1 N N^2], num_elements, num_elements); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.

