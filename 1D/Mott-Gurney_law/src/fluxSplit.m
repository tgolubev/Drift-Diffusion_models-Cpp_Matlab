function [vp,vn] = fluxSplit(u,f,df,strategy)
% WENO flux spliting subroutine.
% OUTPUT:
%   * vp: positive flux, v^{+}, which corresponds to f_{i+1/2}^{-}
%   * vn: negative flux  v^{-}, which corresponds to f_{i+1/2}^{+}

switch strategy
    case{1} % Godunov - scalar fluxsplit (non-conservative)
        vp = f((u + abs(u))./2); %flux^{+}
        vn = f((u - abs(u))./2); %flux^{-}
    case{2} % Local Lax-Friedrichs
        v = f(u); alpha = abs(df(u));
        vp = 0.5.*(v + alpha.*u); %flux^{+}
        vn = 0.5.*(v - alpha.*u); %flux^{-}
    case{3} % (Global) Lax-Friedrichs  %This is most commonly used (see pg. 23 Shu paper: weno for convection dominated problems)
        v = f(u); alpha = max(abs(df(u)));
        vp = 0.5.*(v + alpha.*u); %flux^{+}
        vn = 0.5.*(v - alpha.*u); %flux^{-}
    otherwise
        error('only cases 1,2, and 3 are available')
end