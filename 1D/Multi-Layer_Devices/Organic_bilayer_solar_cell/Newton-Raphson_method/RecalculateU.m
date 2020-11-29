function RecalculateU = RecalculateU(V_left,V_right, n_full, p_full)

global l Vt L_int k_ppr_0 Jx k_ppd_eq manual_Rp E_gap N N_LUMO N_HOMO q Rp_man

F = (V_right-V_left).*Vt./L_int;  %V_left = V(l-1), V_right = V(l) nd recall that V(l) and V(l-1) correspond to fullV(l+1) and fullV(l)
k_ppd =  kppd(F);
Rtotal = R(n_full, p_full);  %calculates recombination:  we use trap and SRH so implement just those 2--> do in  sepearte function file
k_ppr = k_ppr_0*exp(-F.*L_int./Vt); %ISSUE IS HERE!
zeta_eq = Rtotal./k_ppd_eq;   %use ./ to tell it that these are not matrixes! otherwise get error.
zeta = (Jx./L_int+ k_ppr.*zeta_eq + Rtotal)./(k_ppr+k_ppd);
netU = k_ppd*zeta - Rtotal;

if(manual_Rp)
        temp = ( E_gap - F*L_int+ Vt*log(n_full(l+1)*p_full(l)*N^2/(N_LUMO*N_HOMO)) )/(q*L_int*Rp_man);  %IN TEH FIRST TRIAL, N_FULL L+1 IS NOT DEFINED!!! IS JUST = 0!!: NEED TO FIX THIS!!
        netU = netU-temp;
 end

RecalculateU = netU;