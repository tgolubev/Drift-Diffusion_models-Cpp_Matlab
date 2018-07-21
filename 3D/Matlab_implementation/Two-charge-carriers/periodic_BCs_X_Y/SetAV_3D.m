% Setup of 3D Poisson matrix using loop and spdiag

function AV = SetAV_3D(epsilon_avged)

global num_elements N

AV_val = zeros(2*num_elements, 7);   %this is a matrix which will just store the non-zero diagonals of 3D Poisson matrix

eps_X_avg = epsilon_avged.eps_X_avg;
eps_Y_avg = epsilon_avged.eps_Y_avg;
eps_Z_avg = epsilon_avged.eps_Z_avg;



%--------------------------------------------------------------------------
%Lowest diagonal: X's Left PBC
values_cut = zeros(N+1, N+1, N);     %make i and j have N+1 b/c we are including the right bndry (n+1th) pt... for PBCs. Note: k still has N pts, b/c no pbc there
values_cut(1,1:N+1,1:N) = -eps_X_avg(2,2:N+1+1,2:N+1); %just index = 2 for x values--> since all correspond to left X pbc
values_cut = permute(values_cut, [3 2 1]);  
AV_val(1:num_elements, 1) = values_cut(:);  %recall lower diag's fill from the top, so don't need to worry about extra elements etc...


%--------------------------------------------------------------------------

%Lower diagonal:  X's
values_cut = zeros(N+1, N+1, N);     %make i and j have N+1 b/c we are including the right bndry (n+1th) pt... for PBCs. Note: k still has N pts, b/c no pbc there
values_cut(1:N+1-1,1:N+1,1:N) = -eps_X_avg(2+1:N+1+1,2:N+1+1,2:N+1); 
values_cut = permute(values_cut, [3 2 1]);  %this changes order of indices, so have z, y, x --> z are the rows--> then when use (:), will have correct order in AV_val
AV_val(1:num_elements, 2) = values_cut(:); 


%--------------------------------------------------------------------------
%Lower diagonal: Y's LEFT PBC's
values_cut = zeros(N+1, N+1, N);      %RESET THE SIZE OF VALUES_CUT: this is
                                      %important!!                                      
  %note: made the j of N+1 size-->  so elements unfilled --> this corresponds to the 0's subblocks
values_cut(1:N+1,1,1:N) = -eps_Y_avg(2:N+1+1, 2, 2:N+1);  %Note: made y values be just 1 index --> b/c this pbc--> is just always refering to j = 1 left pbc.
values_cut = permute(values_cut, [3 2 1]);
AV_val(1:num_elements, 3) = values_cut(:);  %is N^3 lenght, but the last block will just have 0's....

%--------------------------------------------------------------------------
%Lower diagonal: Y's
values_cut = zeros(N+1, N+1, N);          %RESET THE SIZE OF VALUES_CUT: this is
                                      %important!!                                    
  %note: made the j of N+1 size-->  so when leave the last Nth elements,
  %unfilled--> this corresponds to the 0's subblocks

values_cut(1:N+1,1:N+1-1,1:N) = -eps_Y_avg(2:N+1+1,2+1:N+1+1,2:N+1);  
values_cut = permute(values_cut, [3 2 1]);
AV_val(1:num_elements, 4) = values_cut(:);  %is N^3 lenght, but the last block will just have 0's....

%--------------------------------------------------------------------------
%main lower diagonal (below main diagonal)
values_cut = zeros(N+1, N+1, N); 
values_cut(1:N+1,1:N+1,1:N-1) = -eps_Z_avg(2:N+1+1,2:N+1+1,2+1:N+1);  
values_cut = permute(values_cut, [3 2 1]);
AV_val(1:num_elements, 5) = values_cut(:); 

%--------------------------------------------------------------------------
%main diagonal
values_cut = zeros(N+1, N+1, N); 
values_cut(1:N+1,1:N+1,1:N) = eps_X_avg(2:N+1+1, 2:N+1+1,2:N+1) + eps_X_avg(2+1:N+1+1+1, 2:N+1+1,2:N+1) + eps_Y_avg(2:N+1+1,2:N+1+1,2:N+1) + eps_Y_avg(2:N+1+1,2+1:N+1+1+1,2:N+1) + ...
                              +eps_Z_avg(2:N+1+1,2:N+1+1,2:N+1) + eps_Z_avg(2:N+1+1,2:N+1+1,2+1:N+1+1);  %note: k values go only to N+1 or max to N+1+1...., 1 less than X and Y
                          
values_cut = permute(values_cut, [3 2 1]);
AV_val(1:num_elements,6) = values_cut(:);

%--------------------------------------------------------------------------
%main uppper diagonal
values_cut = zeros(N+1, N+1, N);
values_cut(1:N+1, 1:N+1,1:N-1) = -eps_Z_avg(2:N+1+1, 2:N+1+1, 2+1:N+1);  %2+1 b/c is k+1+1
values_cut = permute(values_cut, [3 2 1]);
AV_val(2:num_elements+1,7) = values_cut(:);  %+1 here b/c of way spdiags fills


%--------------------------------------------------------------------------
%upper diagonal   Y's
values_cut = zeros(N+1, N+1, N);
values_cut(1:N+1, 1:N+1-1,1:N) = -eps_Y_avg(2:N+1+1, 2+1:N+1+1, 2:N+1);  %note is j+1+1
values_cut = permute(values_cut, [3 2 1]);
AV_val(1+N:num_elements+N, 8) = values_cut(:);  %must shift by N+1


%--------------------------------------------------------------------------
%upper diagonal: right Y PBC's
values_cut = zeros(N+1, N+1, N);
values_cut(1:N+1, 1, 1:N) = -eps_Y_avg(2:N+1+1, N+1+1, 2:N+1);  %j = N+1+1 b/c corresponding to right pbc's
values_cut = permute(values_cut, [3 2 1]);
AV_val(1+N^2:num_elements+N^2, 9) = values_cut(:);  %must shift by N+1

%--------------------------------------------------------------------------
%upper diagonal: X's
values_cut = zeros(N+1, N+1, N);
values_cut(1:N+1-1, 1:N+1,1:N) = -eps_X_avg(2+1:N+1+1, 2:N+1+1, 2:N+1); %i+1+1
values_cut = permute(values_cut, [3 2 1]);
AV_val(1+N*(N+1):num_elements + N*(N+1) ,10) = values_cut(:);  %shifted by N*N+1 b/c N z values and N+1 y values to shift over by.

%--------------------------------------------------------------------------
%far upper diagonal: X right PBC's
values_cut = zeros(N+1, N+1, N);
values_cut(1, 1:N+1,1:N) = -eps_X_avg(N+1+1, 2:N+1+1, 2:N+1); %i = N+1+1 b/c always corresonding to far right BC
values_cut = permute(values_cut, [3 2 1]);
AV_val(1+(N^2)*(N+1):num_elements + (N^2)*(N+1) ,11) = values_cut(:);  %shifted b/c spdiags fills upper diags from the bottom

%--------------------------------------------------------------------------
AV = spdiags(AV_val, [-(N^2)*(N+1) -N*(N+1) -N^2 -N -1 0 1 N N^2 N*(N+1) (N^2)*(N+1)], num_elements, num_elements); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.



%=======================OLD VERSION, NO PBC'S==================================================





















% %--------------------------------------------------------------------------
% %Lowest diagonal: corresponds to V(i, j, k-1)
% i = 1;
% j = 1;
% k = 2;
% for index = 1:N^3 - N^2      % (1st element corresponds to Nth row  (number of elements = N^3 - N^2) 
%     AV_val(index,1) = -(epsilon(i+1,j+1, k+1) + epsilon(i+1+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1))/4.;  %spsilons correspond  to full device, so need to +1 from i,j values
%     
%     i = i+1;
%     if (i > N)
%         i = 1;  %reset i when reach end of subblock
%         j = j+1;
%     end
%     if (j > N)
%         j = 1;
%         k = k+1;
%     end
%    
% end
% %--------------------------------------------------------------------------
% %lower diag (y's)
% i = 1;
% j = 2;
% k = 1;
% for index = 1:N^3 - N       
%     if (j > 1)
%         AV_val(index,2) = -(epsilon(i+1,j+1, k+1) + epsilon(i+1+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1))/4.;  %THIS IS WRONG!
%     end
%     
%     i = i+1;
%     if (i > N)
%         i = 1; 
%         j = j+1;
%     end
%     if (j > N)  
%         j = 1;
%         k = k+1;
%     end
%     
%     
% end
% %--------------------------------------------------------------------------
% %main lower diagonal
% i = 2;
% j = 1;
% k = 1;
% for index = 1:N^3 - 1
%     if (i > 1)
%         AV_val(index,3) =  - (epsilon(i+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1, j+1, k+1+1) + epsilon(i+1, j+1+1, k+1+1))/4.;
%     end
%     
%     i = i+1;
%     if (i > N)  
%         i = 1; 
%         j = j+1;
%     end
%     if (j > N)
%         j = 1;
%         k = k+1;
%     end
% end
%     
% %--------------------------------------------------------------------------
% %main diagonal
% i = 1;
% j = 1;
% k = 1;
% for index =  1:num_elements      
%     AV_val(index,4) = ((epsilon(i+1,j+1, k+1) + epsilon(i+1+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1))/4. ...  %... to wrap code
%                     + (epsilon(i+1,j+1, k+1) + epsilon(i+1+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1))/4. ...
%                     + (epsilon(i+1, j+1, k+1) + epsilon(i+1, j+1+1, k+1) + epsilon(i+1, j+1, k+1+1) + epsilon(i+1, j+1+1, k+1+1))/4. ...
%                     + (epsilon(i+1+1,j+1, k+1) + epsilon(i+1+1,j+1+1, k+1) + epsilon(i+1+1, j+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4. ...
%                     + (epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1) + epsilon(i+1, j+1+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4. ...
%                     + (epsilon(i+1,j+1, k+1+1) + epsilon(i+1+1,j+1, k+1+1) +  epsilon(i+1, j+1+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4.);
%                 
%     i = i+1;
%     if (i > N)
%         i = 1;
%         j = j+1;
%     end
%     if (j > N)
%         j = 1;
%         k = k+1;
%     end
% end
% %--------------------------------------------------------------------------
% %main uppper diagonal
% i = 1;
% j = 1;
% k = 1;
% for index = 2:N^3 - 1 +1     % matlab fills this from the bottom (so i = 2 corresponds to 1st row in matrix)
%    if (i > 0)
%         AV_val(index,5) = -(epsilon(i+1+1,j+1, k+1) + epsilon(i+1+1,j+1+1, k+1) + epsilon(i+1+1, j+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4.;
%    end
%    
%    i = i+1;
%    if (i > N-1)  
%        i = 0;
%        j = j+1;
%    end
%     if (j > N)
%         j = 1;
%         k = k+1;
%     end
%    
% end
% %--------------------------------------------------------------------------
% %upper diagonal (y's)
% i = 1;
% j = 1;
% k = 1;
% for index = 1+N:N^3 - N +N  %matlab fills this from bottom... (in cpp will just start from 1)
%     if (j > 0)
%         AV_val(index, 6) = -(epsilon(i+1, j+1+1, k+1) + epsilon(i+1+1, j+1+1, k+1) + epsilon(i+1, j+1+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4.;
%     end
%     
%     i = i+1;
%      if (i > N)
%         i = 1;
%         j = j+1;
%      end
%      if (j > N-1)
%         j = 0;
%         k = k+1;
%      end
% end
%     
% %--------------------------------------------------------------------------
% %far upper diagonal
% i = 1;
% j = 1;
% k = 1;
% for index = 1+N*N:num_elements  %note: N^3-N^2 elements, but shifted by N^2 b/c of matlab filling convention     %matlab fills from bottom, so this starts at 1+N (since 1st element is in the 2nd subblock of matrix)
%     AV_val(index,7) = -(epsilon(i+1,j+1, k+1+1) + epsilon(i+1+1,j+1, k+1+1) +  epsilon(i+1, j+1+1, k+1+1) + epsilon(i+1+1, j+1+1, k+1+1))/4.;            %1st element corresponds to 1st row.   this has N^2 -N elements
%     
%     i = i+1;
%     if (i > N)
%         i = 1;
%         j = j+1;
%     end
%     if (j > N)
%         j = 1;
%         k = k+1;
%     end
%     
% end
% 
% %all not specified elements will remain zero, as they were initialized
% %above.
% 
% AV = spdiags(AV_val, [-N^2 -N -1 0 1 N N^2], num_elements, num_elements); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
