function G = GenerationRate()

global G_max num_elements

% G = ones(N,N);
G = ones(num_elements,1);          %better to use it as a 1D vector

%Here can use any relevant analytic expression or input data from a file
% i.e.
%G(2:num_cell) = load('gen_rate.txt');
%G(x,z)
% G(1:N,1:N) = 1;  %these just include the insides--> since want Un and Up to just include the  insides--> i,j directly correspond to the matrix indices

G = G_max.*G;  %skip the endpoints b/c we don't need them in Un...




