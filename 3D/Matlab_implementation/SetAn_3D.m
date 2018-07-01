% Setup of 3D electron continuity matrix using loop and spdiag
%Note: this is almost the same as Ap setup, except the Bn_pos and Bn_neg
% are switched everywhere

function An = SetAn_3D(n_mob_avged, Bernoulli_n_values)

global num_elements N

%extract variables from struct, for brevity in eqns
Bn_posX = Bernoulli_n_values.Bn_posX;
Bn_negX = Bernoulli_n_values.Bn_negX;
Bn_posY = Bernoulli_n_values.Bn_posY;
Bn_negY = Bernoulli_n_values.Bn_negY;
Bn_posZ = Bernoulli_n_values.Bn_posZ;
Bn_negZ = Bernoulli_n_values.Bn_negZ;

n_mob_X_avg = n_mob_avged.n_mob_X_avg;
n_mob_Y_avg = n_mob_avged.n_mob_Y_avg;
n_mob_Z_avg = n_mob_avged.n_mob_Z_avg;

An_val = zeros(num_elements, 7);   %this is a matrix which will just store the non-zero diagonals of 2D hole continuity eqn

%NOTE: index is not neccesarily equal to i (the x index of V), it is the
%index of the diagonal arrays.
%--------------------------------------------------------------------------
%Lowest diagonal: 
values = -n_mob_Z_avg.*Bn_negZ;
values_cut = zeros(N, N, N);
values_cut(1:N,1:N,1:N-1) = values(2:N+1,2:N+1,2+1:N+1);  %shift by +1,   %note: the k values  uses 2+1:N+1 -->  b/c k = 2:N is what is needed..., so need to skip k = 1, which corresponds to 2 in orig values matrix
An_val(1:N^3, 1) = values_cut(:);

% index = 1;
% for k = 2:N
%     for j = 1:N
%         for i = 1:N
%              An_val(index,1) = values(i+1,j+1,k+1);
%             index = index +1;
%         end 
%     end
% end

%--------------------------------------------------------------------------
%lower diagonal
values = -n_mob_Y_avg.*Bn_negY;
values_cut = zeros(N, N, N);          %RESET THE SIZE OF VALUES_CUT: this is
                                      %important!!
                                      
  %note: made the j of N size-->  so when leave the last Nth elements,
  %unfilled--> this corresponds to the 0's subblocks

values_cut(1:N,1:N-1,1:N) = values(2:N+1,2+1:N+1,2:N+1);  
An_val(1:N^3, 2) = values_cut(:); 

% index = 1;
% for k = 1:N
%     for j = 2:N
%         for i = 1:N
%             An_val(index,2) = values(i+1,j+1,k+1);
%             index = index+1;
%         end
%     end
%     index = index + N;  %add on the 0's subblock
% end

%--------------------------------------------------------------------------
%main lower diagonal (below main diagonal)
values = -n_mob_X_avg.*Bn_negX;
values_cut = zeros(N, N, N); %again made i have N elements here, but fill to N-1--> to have the  0's in corners where need them
values_cut(1:N-1,1:N,1:N) = values(2+1:N+1,2:N+1,2:N+1); 
An_val(1:N^3, 3) = values_cut(:);


% index = 1;
% for k = 1:N
%     for j = 1:N
%         for i = 2:N
%             An_val(index,3) = values(i+1,j+1,k+1);
%             index = index+1;
%         end
%         index = index +1;  %add on the corner element which is 0
%     end
% end
   
%--------------------------------------------------------------------------
%main diagonal
values1 = n_mob_Z_avg.*Bn_posZ + n_mob_Y_avg.*Bn_posY + n_mob_X_avg.*Bn_posX;
values2 = n_mob_X_avg.*Bn_negX;
values3 = n_mob_Y_avg.*Bn_negY;
values4 = n_mob_Z_avg.*Bn_negZ;

% values_cut = zeros(N,N,N);  %preallocation is not neccessary  here, since
% all the values are filled
values_cut(1:N,1:N,1:N) = values1(2:N+1, 2:N+1,2:N+1) + values2(2+1:N+1+1, 2:N+1,2:N+1) + values3(2:N+1,2+1:N+1+1,2:N+1) + values4(2:N+1,2:N+1,2+1:N+1+1);
An_val(1:N^3,4) = values_cut(:);

% index = 1;
% for k = 1:N
%     for j = 1:N
%         for i = 1:N     
%             An_val(index,4) = values1(i+1,j+1,k+1) + values2(i+1+1,j+1,k+1) + values3(i+1,j+1+1,k+1) + values4(i+1,j+1,k+1+1);
%             index = index+1;
%         end
%     end
% end

%--------------------------------------------------------------------------
%main uppper diagonal
values = -n_mob_X_avg.*Bn_posX;
values_cut = zeros(N,N,N);
values_cut(1:N-1, 1:N,1:N) = values(2:N, 2:N+1, 2:N+1);

An_val(2:N^3+1,5) = values_cut(:);



% index = 2;
% for k = 1:N
%     for j = 1:N
%         for i = 1:N-1     
%             An_val(index,5) = values(i+1+1,j+1,k+1);
%             index = index+1;
%         end
%         index = index+1;  %add on the element which = 0 --> corner...
%     end
% end
       
%--------------------------------------------------------------------------
%upper diagonal
values = -n_mob_Y_avg.*Bn_posY;
values_cut = zeros(N,N,N);
values_cut(1:N, 1:N-1,1:N) = values(2:N+1, 2:N, 2:N+1);

An_val(1+N:N^3+N, 6) = values_cut(:);

% index = 1+N;
% for k = 1:N
%     for j = 1:N-1
%         for i = 1:N    
%             An_val(index, 6) = values(i+1,j+1+1,k+1);
%             index = index+1;
%         end
%     end
%     index = index + N;  %add on the empty block--> of 0's corresponding to y BC's. 
% end
            
%--------------------------------------------------------------------------
%far upper diagonal
values = -n_mob_Z_avg.*Bn_posZ;
values_cut = zeros(N,N,N);
values_cut(1:N, 1:N,1:N-1) = values(2:N+1, 2:N+1, 2:N);

An_val(1+N*N:N^3+N^2,7) = values_cut(:);

% index = 1+N*N;
% for k = 1:N-1
%     for j = 1:N
%         for i = 1:N    
%              An_val(index,7) = values(i+1,j+1,k+1+1);  
%              index = index+1;
%         end
%     end
% end


%all not specified elements will remain zero, as they were initialized
%above.


An = spdiags(An_val, [-N^2 -N -1 0 1 N N^2], num_elements, num_elements); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.

