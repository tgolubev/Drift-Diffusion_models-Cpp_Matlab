% Setup of 3D Poisson matrix using loop and spdiag

function AV = SetAV_3D(epsilon)

global num_elements N

AV_val = zeros(num_elements, 7);   %this is a matrix which will just store the non-zero diagonals of 3D Poisson matrix

%NOTE: index is not neccesarily equal to i (the x index of V), it is the
%index of the diagonal arrays.
%--------------------------------------------------------------------------
%Lowest diagonal: corresponds to V(i, j, k-1)
i = 1;
j = 1;
k = 2;
for index = 1:N^3 - N^2      % (1st element corresponds to Nth row  (number of elements = N^3 - N^2) 
    AV_val(index,1) = -(epsilon(i+1,j+1, k+1) + epsilon(i+1+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1))/4.;  %spsilons correspond  to full device, so need to +1 from i,j values
    
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
%lower diag (y's)
i = 1;
j = 2;
k = 1;
for index = 1:N^3 - N       
    if (j > 1)
        AV_val(index,2) = -(epsilon(i+1,j+1, k+1) + epsilon(i+1+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1))/4.;
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
%main lower diagonal
i = 2;
j = 1;
k = 1;
for index = 1:N^3 - 1
    if (i > 1)
        AV_val(index,3) =  - (epsilon(i+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1, j+1, k+1+1) + epsilon(i+1, j+1+1, k+1+1))/4.;
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
    AV_val(index,4) = ((epsilon(i+1,j+1, k+1) + epsilon(i+1+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1))/4. ...  %... to wrap code
                    + (epsilon(i+1,j+1, k+1) + epsilon(i+1+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1))/4. ...
                    + (epsilon(i+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1, j+1, k+1+1) + epsilon(i+1, j+1+1, k+1+1))/4. ...
                    + (epsilon(i+1+1,j+1, k+1) + epsilon(i+1+1,j+1+1, k+1) + epsilon(i+1+1, j+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4. ...
                    + (epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1) + epsilon(i+1, j+1+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4. ...
                    + (epsilon(i+1,j+1, k+1+1) + epsilon(i+1+1,j+1, k+1+1) +  epsilon(i+1, j+1+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4.);
                
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
for index = 2:N^3 - 1 +1     % matlab fills this from the bottom (so i = 2 corresponds to 1st row in matrix)
   if (i > 0)
        AV_val(index,5) = -(epsilon(i+1+1,j+1, k+1) + epsilon(i+1+1,j+1+1, k+1) + epsilon(i+1+1, j+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4.;
   end
   
   i = i+1;
   if (i > N-1)  
       i = 0;
       j = j+1;
   end
    if (j > N)
        j = 1;
        k = k+1;
    end
   
end
%--------------------------------------------------------------------------
%upper diagonal (y's)
i = 1;
j = 1;
k = 1;
for index = 1+N:N^3 - N +N  %matlab fills this from bottom... (in cpp will just start from 1)
    if (j > 0)
        AV_val(index, 6) = -(epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1) + epsilon(i+1, j+1+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4.;
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
for index = 1+N*N:num_elements  %note: N^3-N^2 elements, but shifted by N^2 b/c of matlab filling convention     %matlab fills from bottom, so this starts at 1+N (since 1st element is in the 2nd subblock of matrix)
    AV_val(index,7) = -(epsilon(i+1,j+1, k+1+1) + epsilon(i+1+1,j+1, k+1+1) +  epsilon(i+1, j+1+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4.;            %1st element corresponds to 1st row.   this has N^2 -N elements
    
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

AV = spdiags(AV_val, [-N^2 -N -1 0 1 N N^2], num_elements, num_elements); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
