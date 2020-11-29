function Fn_element = Continuity_n_fnc(Jn_r, Jn_l, Uni)
%Jn_r and Jn_l stand for right and left
global dx q Cn

Fn_element = (Jn_r - Jn_l) + Cn*Uni;              %NOTE: Jn is defined on right side: i.e. Jn_{i-1/2} is Jn_i
