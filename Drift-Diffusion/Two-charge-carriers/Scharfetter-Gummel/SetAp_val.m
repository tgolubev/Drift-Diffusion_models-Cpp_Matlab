function Ap_val = SetAp_val(num_elements, B, fullV, p, Vt) 

Ap =  zeros(num_elements,3);

                
  for i=1:num_elements
    Ap(i,1) = B(1,i);                        %lower diagonal
    Ap(i,2) = -(B(2,i) + B(1,i+1));          %main diagonal
    Ap(i,3) = B(2,i+1);                      %upper diagonal   CHECK MAKE SURE THIS IS RIGHT!
 end

 
Ap_val = spdiags(Ap,-1:1,num_elements,num_elements); %A = spdiags(B,d,m,n) creates an m-by-