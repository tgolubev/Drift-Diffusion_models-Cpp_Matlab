function G = GenerationRate()

global G_max num_cell

%Here can use any relevant analytic expression or input data from a file
% i.e.
%G(2:num_cell) = load('gen_rate.txt');

G(1:num_cell-1, 1:num_cell-1) = 1;  %ONLY DEFINE IT IN THE INTERIOR OF THE DEVICE

G = G_max*G;  %skip the endpoints b/c we don't need them in Un...




