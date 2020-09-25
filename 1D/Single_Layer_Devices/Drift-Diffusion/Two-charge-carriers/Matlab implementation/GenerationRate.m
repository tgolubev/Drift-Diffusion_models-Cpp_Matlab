function G = GenerationRate()

global G_max num_cell

%Here can use any relevant analytic expression or input data from a file
% i.e.
%G = load('gen_rate.inp'); %skip the endpoints b/c we don't need them in Un....

 G(2:num_cell) = 1;

 G = G_max*G(2:num_cell)/max(G(2:num_cell));  %normalize it  




