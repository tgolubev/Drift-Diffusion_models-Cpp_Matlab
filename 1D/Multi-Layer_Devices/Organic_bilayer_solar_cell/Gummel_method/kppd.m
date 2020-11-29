function k_ppd = kppd(F)
global q epsilon Vt Eb k_rec L_int 


if(F <= 0.0D0)    %This performs the integral over half sphere for case when field opposes dissociation
    zeta1 = -pi/2.0D0;
    zeta2 = pi/2.0D0;
    n = 1000;
    sumk = 0.0D0;
    for k=1:n-1
        sumk = sumk + f_of_zeta(zeta1 + k*(zeta2-zeta1)/n, F);
    end 
    k_ppd = (zeta2-zeta1)./n.*(f_of_zeta(zeta1, F)./2.0D0 + sumk + f_of_zeta(zeta1, F)./2.0D0);
    k_ppd = k_ppd./pi;
else
 
    bk = q*F/(8*pi*epsilon*(Vt^2));  %b in Bessel fnc. experssion (**2 means ^2)   %USE ./ to tell matlab that its element by element operation --> BUT THIS THEN RETURNS AN EMPTY MATRIX!
    k_ppd = (3.0.*k_rec./(4.0*pi.*L_int^3))*exp(-Eb./Vt)*(1.0 + bk + (bk^2)/3.0 + (bk^3)/18.0 + (bk^4)/180.0);  %Bessel function is expanded to 4th order here
end

%not sure what these are for
%bk = q*F/(8.0*pi*epsilon*Vt^2)
%k_ppd2 = (3.0*k_rec/(4.0*pi*L_int^3))*exp(-Eb/Vt)*(1.0 + bk + (bk^2)/3.0 + (bk^3)/18.0D0 + (bk^4)/180.0D0);  %Bessel function is expanded to 4th order here