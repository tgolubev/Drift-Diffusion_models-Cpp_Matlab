function G = GenerationRate(gen_rate_input)
global l_HTL_int l_ETL_int num_cell G_max

% This uses real generation rate from optical model. 
% Reads in data from file gen_rate.txt (must be in same directory)
% gen_rate.txt contains generation rate at each mesh grid point in the
% simulation. The # of points in the file must match the mesh
% used in the code.

G(1:l_HTL_int) = 0;   % no charge generation in HTL

G(l_HTL_int+1:l_ETL_int) = gen_rate_input;

%normalize by the max value
G(l_HTL_int+1:l_ETL_int) =  G(l_HTL_int+1:l_ETL_int)./max(G(l_HTL_int+1:l_ETL_int));

G(l_ETL_int+1:num_cell+1) = 0;  % no charge generation in ETL

G = G_max*G(2:num_cell);  %skip the endpoints b/c we don't need thme in Un....



