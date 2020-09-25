function Ap_val = SetAp_val(B_p) 
global num_elements p_mob 

Ap =  zeros(num_elements,3);

for i=1:num_elements     
    Ap(i,2) = -(p_mob(i+1)*B_p(2,i+1) + p_mob(i+1+1)*B_p(1,i+1+1));    
end

for i = 1:num_elements-1
    Ap(i,1) = p_mob(i+1+1)*B_p(1,i+1+1);   
end
Ap(num_elements, 1) = 0;     % last element is unused
for i = 2:num_elements
    Ap(i,3) = p_mob(i+1)*B_p(2,i+1);     
end

Ap(1,3) = 0;                 % 1st element of Ap(:,3) is unused

Ap_val = spdiags(Ap,-1:1,num_elements,num_elements);
