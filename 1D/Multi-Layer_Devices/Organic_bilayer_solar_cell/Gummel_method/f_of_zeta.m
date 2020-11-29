function fnc_of_zeta = f_of_zeta(zeta, F)  %make returned value have diff name than fnc, to not have confusion
global L_int k_rec Eb rc Vt 
fnc_of_zeta = 3.0D0./(4.0D0.*pi.*L_int^3).*k_rec.*exp((-Eb + F.*rc.*cos(zeta))./Vt); %the integrand