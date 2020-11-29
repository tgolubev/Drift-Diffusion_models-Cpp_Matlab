function Jac = FindJacobian(V, n, p, fullV, n_full, p_full, Un, Up, B_n, B_p)  %note Jacobian is a reserved fnc in matlab
global num_elements N N_LUMO N_HOMO T_tA T_tD CV CnT CpT T dx num_cell l Cp Cn 

%NOTE: Fv, Fn, Fp, Fu will be the 4 equations which will take derivatives
%of

Jac = zeros(3*num_elements, 3*num_elements);


%% F1 = Fv derivatives (Poisson eqn.)

%fill 1st row
Jac(1,1)  = -2.;
Jac(1,2) = 1.;
%fill center part in
 for i = 2:num_elements -1 %i labels matrix columns  %NOTE: HERE, might be ok, if have a lot of 0's --> b/c if V is linear, then the gradient of V is 0.
      Jac(i,i) = -2.;
      Jac(i,i+1) = 1.;
      Jac(i,i-1) = 1.;
 end
 %fill last row
 Jac(num_elements, num_elements-1) = 1.;
 Jac(num_elements, num_elements) = -2.;
 
 
 %NOTE: the 0's DO NOT need to be filled in, b/c i  initiallize the
 %jacobian with  zeros...
 
 %dFv/dn  --> all analytic derivatives
 for i = num_elements +1:2*num_elements  %WORKS CORRECTLY
    Jac(i-num_elements, i) =  -CV;  %diagonal terms are just constant  %i - num_elements labels the columns --> should still be in the 1st section of matrix   
%     Jac(i-num_elements,i)
 end
 %off diag. terms don't nee to be filled b/c are just 0
        
 %dFv/dp --> all analytic
 for i = 2*num_elements +1:3*num_elements
    Jac(i-2*num_elements, i) =  CV;  %diagonal terms are just constant  %i - num_elements labels the columns --> should still be in the 1st section of matrix   
 end

 
 %Fv doesn't depend on U so no extra derivatives needed when add U coupling
 %to Jacobian
 
 
%  %dFv/dnT --> all analytic
%   for i = 3*num_elements +1:4*num_elements
%     Jac(i-3*num_elements, i) =  -CV;  %diagonal terms are just constant  %i - num_elements labels the columns --> should still be in the 1st section of matrix   
%   end
%  
%   %dFv/dpT --> all analytic
%   for i = 4*num_elements +1:5*num_elements
%     Jac(i-4*num_elements, i) =  CV;  %diagonal terms are just constant  %i - num_elements labels the columns --> should still be in the 1st section of matrix   
%  end
 
 %% F2 = Fn derivatives (Jn continuity equation)
 %Here will use just a constant eps_n --> since using different eps_n's for
 %each point is too expensive--> need to recalculate Bernoulli's each
 %time..
 
 %dFn_dVs
 %compute necessary bernoulli derivs w.r.t dVs
  dBni_dVis =dBni_dVi(fullV);
  dnegBni_dVis =  dnegBni_dVi(fullV);
  dBni_dVi_min1s = dBni_dVi_min1(fullV);
  dnegBni_dVi_min1s = dnegBni_dVi_min1(fullV);
 
  dBni_plus1_dVi_plus1(3:num_cell+1) = dBni_dVis(3:num_cell +1);  %is just the dB(dVi)/dVi derivative but replace dVi by dV_i+1d
  dnegBni_plus1_dVi_plus1(3:num_cell+1) = dnegBni_dVis(3:num_cell +1);  %is just the dB(dVi)/dVi derivative but replace dVi by dV_i+1
  dBni_plus1_dVi(3:num_cell+1) =  dBni_dVi_min1s(3:num_cell+1);  %is just a shift of other derivative..
  dnegBni_plus1_dVi(3:num_cell+1) = dnegBni_dVi_min1s(3:num_cell +1);  
   

Jac(1+num_elements,1) = n_full(1)*dnegBni_dVis(2) - n_full(2)*(dBni_dVis(2)  +  dnegBni_plus1_dVi(3)) + n_full(3)*dBni_plus1_dVi(3); %do in terms of n_full b/c need boundary value..
Jac(1+num_elements, 2) = -n(1)*dnegBni_plus1_dVi_plus1(2) + n(2)*dBni_plus1_dVi_plus1(2);
%fill middle
for i = 2:num_elements-1  
       Jac(i+num_elements,i-1) = n(i-1)*dnegBni_dVi_min1s(i) - n(i)*dBni_dVi_min1s(i);
       Jac(i+num_elements,i) =  n(i-1)*dnegBni_dVis(i) - n(i)*(dBni_dVis(i)  +  dnegBni_plus1_dVi(i+1)) + n(i+1)*dBni_plus1_dVi(i+1);
       Jac(i+num_elements,i+1) = -n(i)*dnegBni_plus1_dVi_plus1(i+1) + n(i+1)*dBni_plus1_dVi_plus1(i+1);
 end
 %fill last row
 Jac(2*num_elements,num_elements-1) = n(num_elements-1)*dnegBni_dVi_min1s(num_elements) - n(num_elements)*dBni_dVi_min1s(num_elements);
 Jac(2*num_elements, num_elements) = n_full(num_cell-1)*dnegBni_dVis(num_cell) - n_full(num_cell)*(dBni_dVis(num_cell)  +  dnegBni_plus1_dVi(num_cell+1)) + n_full(num_cell+1)*dBni_plus1_dVi(num_cell+1);

 
 %Modify the  elements which have non-zero Un in them to reflect the
 %partial deriv's of Un....
 epsV = 0.1*abs(V(l) - V(l-1));
 Jac(num_elements + l,num_elements + l-1) =  Jac(num_elements + l,num_elements + l-1) + Cn*(RecalculateU(V(l-1)+epsV/2.,V(l), n_full, p_full) - RecalculateU(V(l-1)-epsV/2., V(l), n_full, p_full))/epsV;
 Jac(num_elements + l,num_elements + l) =   Jac(num_elements + l,num_elements + l) + Cn*(RecalculateU(V(l-1), V(l)+epsV/2., n_full, p_full) - RecalculateU(V(l-1), V(l)-epsV/2., n_full, p_full))/epsV;
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%

%dFn/dns
%first row
Jac(num_elements +1, num_elements +1) = -(B_n(1,2) + B_n(2,3));
Jac(num_elements +1, num_elements +2) = B_n(1,3);
 for i = num_elements +2:2*num_elements-1
     Jac(i,i-1) = B_n(2,i-num_elements+1);
     Jac(i, i) = -(B_n(1,i-num_elements+1) + B_n(2,i+1-num_elements +1));  %here using i,i, since corresponds to block: 2,2 in matrix.  And use i-num_elements, since derivative array has size num_elements...
     Jac(i,i+1) = B_n(1,i+1-num_elements +1);
 end  %using i-num_elements+1  b/c need to subtract off num_elements to be using real device node indices, but B's defined starting from 2.
 %last row
 Jac(2*num_elements, 2*num_elements -1) = B_n(2,num_elements+1);
 Jac(2*num_elements, 2*num_elements) = B_n(1,num_elements +2);
 %recall B_n or B_p(1,i) = B(+dV) and B_p(2,i) = B(-dV)
 %-----------------------------
 %modify elements which have Un in them
 n_perturb = zeros(1, num_elements +2);  %size of n_full
 %so d/dn_l-1 will have n_full(l) modified:
  eps_n = 0.1*abs(n(l) - n(l-1));  %just pick difference between 2 points to get right order of maginitude...
  n_perturb(l) = eps_n/2.;  %NOTE; perturb has size of n_full so use n_full indices
  Jac(num_elements +l, num_elements + l-1) =  Jac(num_elements +l, num_elements + l-1) + Cn*(RecalculateU(V(l-1), V(l), n_full + n_perturb, p_full) - RecalculateU(V(l-1), V(l), n_full - n_perturb, p_full))/eps_n;
  %reset perturbation
  n_perturb = zeros(1, num_elements +2);  %size of n_full
  n_perturb(l+1) = eps_n/2.;  %NOTE; perturb has size of n_full so use n_full indices
  Jac(num_elements +l, num_elements + l) =  Jac(num_elements +l, num_elements + l) + Cn*(RecalculateU(V(l-1), V(l), n_full + n_perturb, p_full) - RecalculateU(V(l-1), V(l), n_full - n_perturb, p_full))/eps_n;
 
 
  %now dFn/dp --> are all 0 EXCEPT for the 2 terms corresponding where have
  %a Un --> Un depends on pfull at l and l+1:
  p_perturb = zeros(1, num_elements +2);  %size of p_full
  %so d/dp_l-1 will have p_full(l) modified:
  eps_p = 0.1*abs(p(l) - p(l-1));  %just pick difference between 2 points to get right order of maginitude...
  p_perturb(l) = eps_p/2.;
  Jac(num_elements + l, 2*num_elements + l-1) = Cn*(RecalculateU(V(l-1), V(l), n_full, p_full+ p_perturb) - RecalculateU(V(l-1), V(l), n_full, p_full-p_perturb))/eps_p;
  p_perturb(l+1) = eps_p/2.;
   Jac(num_elements + l, 2*num_elements + l) = Cn*(RecalculateU(V(l-1), V(l), n_full, p_full+ p_perturb) - RecalculateU(V(l-1), V(l), n_full, p_full-p_perturb))/eps_p;
  

 %all  dFn/dnT, and dFn/dpT are zero --> nothing needs to be filled in

 
 %DONE
 
 %% F3 = Fp derivatives (Jp continuity equation)
 %Here will use just a constant eps_p --> since using different eps_p's for
 %each point is too expensive--> need to recalculate Bernoulli's each
 %time..
 
  dBpi_dVis =dBpi_dVi(fullV);
  dnegBpi_dVis =  dnegBpi_dVi(fullV);
  dBpi_dVi_min1s = dBpi_dVi_min1(fullV);
  dnegBpi_dVi_min1s = dnegBpi_dVi_min1(fullV);
    
  dBpi_plus1_dVi_plus1(3:num_cell+1) = dBpi_dVis(3:num_cell +1);  %is just the dB(dVi)/dVi derivative but replace dVi by dV_i+1
  dnegBpi_plus1_dVi_plus1(3:num_cell+1) = dnegBpi_dVis(3:num_cell +1);  %is just the dB(dVi)/dVi derivative but replace dVi by dV_i+1
  dBpi_plus1_dVi(3:num_cell+1) =  dBpi_dVi_min1s(3:num_cell+1);  %is just a shift of other derivative..
  dnegBpi_plus1_dVi(3:num_cell+1) = dnegBpi_dVi_min1s(3:num_cell +1);  
  
  
Jac(1+2*num_elements,1) = -(p_full(1)*dBpi_dVis(2) - p_full(2)*(dnegBpi_dVis(2)  +  dBpi_plus1_dVi(3)) + p_full(3)*dnegBpi_plus1_dVi(3)); %do in terms of n_full b/c need boundary value..
Jac(1+2*num_elements, 2) = -(-p(1)*dBpi_plus1_dVi_plus1(2) + p(2)*dnegBpi_plus1_dVi_plus1(2));
%fill middle
for i = 2:num_elements-1  
       Jac(i+2*num_elements,i-1) = -(p(i-1)*dBpi_dVi_min1s(i) - p(i)*dnegBpi_dVi_min1s(i));
       Jac(i+2*num_elements,i) =  -(p(i-1)*dBpi_dVis(i) - p(i)*(dnegBni_dVis(i)  +  dBpi_plus1_dVi(i+1)) + p(i+1)*dnegBpi_plus1_dVi(i+1));
       Jac(i+2*num_elements,i+1) = -(-p(i)*dBpi_plus1_dVi_plus1(i+1) + p(i+1)*dnegBpi_plus1_dVi_plus1(i+1));
 end
 %fill last row
 Jac(3*num_elements,num_elements-1) = -(p(num_elements-1)*dnegBpi_dVi_min1s(num_elements) - p(num_elements)*dBpi_dVi_min1s(num_elements));
 Jac(3*num_elements, num_elements) = -(p_full(num_cell-1)*dBpi_dVis(num_cell) - p_full(num_cell)*(dnegBpi_dVis(num_cell)  +  dBpi_plus1_dVi(num_cell+1)) + p_full(num_cell+1)*dnegBpi_plus1_dVi(num_cell+1));

 %modify elements corresponding to Up(l-1)--> for d/dV...
  epsV = 0.1*abs(V(l) - V(l-1));
  Jac(2*num_elements + l-1,num_elements + l-1) =  Jac(2*num_elements + l-1,num_elements + l-1) - Cp*(RecalculateU(V(l-1)+epsV/2.,V(l), n_full, p_full) - RecalculateU(V(l-1)-epsV/2., V(l), n_full, p_full))/epsV;
  Jac(2*num_elements + l-1,num_elements + l) =   Jac(2*num_elements + l-1,num_elements + l) - Cp*(RecalculateU(V(l-1), V(l)+epsV/2., n_full, p_full) - RecalculateU(V(l-1), V(l)-epsV/2., n_full, p_full))/epsV;
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%recall that Bernoulli fnc is independent of n and p, so here do NOT need
 %to recalculate Bernoulli's..., but do need to recalculate the J's
%  dFp_dps = dFp_dp(B_p, p_full, Up, eps_p);
%first row
Jac(2*num_elements +1, 2*num_elements +1) = (B_p(2,2) + B_n(1,3));
Jac(2*num_elements +1, 2*num_elements +2) = -B_p(2,3);
 for i = 2*num_elements +2:3*num_elements-1
     Jac(i,i-1) = -B_p(1,i-2*num_elements+1);
     Jac(i, i) = (B_p(2,i-2*num_elements+1) + B_p(1,i+1-2*num_elements +1));  %here using i,i, since corresponds to block: 2,2 in matrix.  And use i-num_elements, since derivative array has size num_elements...
     Jac(i,i+1) = -B_p(2,i+1-2*num_elements +1);
 end  %using i-num_elements+1  b/c need to subtract off num_elements to be using real device node indices, but B's defined starting from 2.
 %last row
 Jac(3*num_elements, 3*num_elements -1) = -B_p(1,num_elements+1); %NOTE: using num_elements here b/c B_p's etc defined only on 1 device lenght indices....
 Jac(3*num_elements, 3*num_elements) = (B_p(2,num_elements +1) + B_p(1,num_elements +2));
 
 %modify elements which have Up in them: Up is at l-1: so dFp{l-1} terms
 %need to be modified
 n_perturb = zeros(1, num_elements +2);  %size of n_full
 %so dFp/dn ARE ALL 0 EXCEPT where have Up --> at l-1
  eps_n = 0.1*abs(n(l) - n(l-1));  %just pick difference between 2 points to get right order of maginitude...
  n_perturb(l) = eps_n/2.;  %NOTE; perturb has size of n_full so use n_full indices
  Jac(2*num_elements +l-1, num_elements + l-1) =  - Cp*(RecalculateU(V(l-1), V(l), n_full + n_perturb, p_full) - RecalculateU(V(l-1), V(l), n_full - n_perturb, p_full))/eps_n;
  %reset perturbation
  n_perturb = zeros(1, num_elements +2);  %size of n_full
  n_perturb(l+1) = eps_n/2.;  %NOTE; perturb has size of n_full so use n_full indices
  Jac(2*num_elements +l-1, num_elements + l) =  - Cp*(RecalculateU(V(l-1), V(l), n_full + n_perturb, p_full) - RecalculateU(V(l-1), V(l), n_full - n_perturb, p_full))/eps_n;
 
 
  %now dFp/dp --> 
  %a Un --> Un depends on pfull at l and l+1:
  p_perturb = zeros(1, num_elements +2);  %size of p_full
  %so d/dp_l-1 will have p_full(l) modified:
  eps_p = 0.1*abs(p(l) - p(l-1));  %just pick difference between 2 points to get right order of maginitude...
  p_perturb(l) = eps_p/2.;
  Jac(2*num_elements + l-1, 2*num_elements + l-1) =  Jac(2*num_elements + l-1, 2*num_elements + l-1) -  Cp*(RecalculateU(V(l-1), V(l), n_full, p_full+ p_perturb) - RecalculateU(V(l-1), V(l), n_full, p_full-p_perturb))/eps_p;
  p_perturb(l+1) = eps_p/2.;
   Jac(2*num_elements + l-1, 2*num_elements + l) =  Jac(2*num_elements + l-1, 2*num_elements + l) - Cp*(RecalculateU(V(l-1), V(l), n_full, p_full+ p_perturb) - RecalculateU(V(l-1), V(l), n_full, p_full-p_perturb))/eps_p;
  
 
 


 %all dFp/dn, dFp/dnT, and dFp/dpT are zero --> nothing needs to be filled in
 
 
 %NOTE: THERE ARE NO TRAPS CONSIDERED IN THIS VERSION!! SO NO F4 and F5
 %functions
 
%  %% F4 = FnT derivatives (nT trap eqn)
%  %all derivatives are analytic here
%  %Only blocks 4,2 and 4,4 have non-zero terms
%  
%  %block 4,2 (row, column)
%  for i = num_elements +1: 2*num_elements
%      Jac(i+2*num_elements,i) = -CnT*(n(i-num_elements)*N/N_LUMO)^(T/T_tA-1);
%  end
%  
%  %block 4,4
%  for i = 3*num_elements +1: 4*num_elements
%      Jac(i,i) = 1.0;  %here; derivative of FnT (= nT) is just 1
%  end
%  
% % dFnT/dV, dFnT/dp, and dFnT/dpT are all zero --> nothing needs to be
% % filled in
%  
%  %% F5 = FpT derivatives (pT trap eqn)
%  %all derivatives are analytic here
%  %Only blocks 5,3 and 5,5 have non-zero terms
%  
%  %block 5,3 (row, column)
%  for i = 2*num_elements +1: 3*num_elements
%      Jac(i+2*num_elements,i) = -CpT*(p(i-2*num_elements)*N/N_HOMO)^(T/T_tD-1);
%  end
%  
%  %block 5,5
%  for i = 4*num_elements +1: 5*num_elements
%      Jac(i,i) = 1.0;  %here; derivative of FnT (= nT) is just 1
%  end
 
% dFpT/dV, dFpT/dn, and dFpT/dnT are all zero --> nothing needs to be
% filled in
