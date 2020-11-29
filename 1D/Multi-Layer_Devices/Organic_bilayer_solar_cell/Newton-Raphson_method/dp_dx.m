function dp_dx = dp_dx(p_full)
global dx num_elements

for i =  2:num_elements+1  %shift by 1, b/c boundary values of n_full start at 1
    dp_dx(i-1) = (p_full(i+1) - p_full(i-1))/(2*dx);  %shift to i-1 b/c dn_dx defines inside values only
end
