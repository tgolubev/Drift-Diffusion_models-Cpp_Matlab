function Fp_element = Continuity_p_fnc(Jp_r, Jp_l, Upi)  

global dx q Cp

Fp_element = (Jp_r - Jp_l) - Cp*Upi;              %NOTE: Jn is defined on right side: i.e. Jn_{i-1/2} is Jn_i
