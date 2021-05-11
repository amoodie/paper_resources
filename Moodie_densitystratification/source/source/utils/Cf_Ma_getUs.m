% Based on Ma et at. Resistance relation...
% Calculate shear velocity Us given bulk velocity and bed grain size
function [Us] = Cf_Ma_getUS(U,vs)

    N=length(U);         % replace it with number of data. e.g. length(U)
    
    Us0=0.05*ones(N,1);  % Initial values, should not matter.
    
    Us = fsolve(@(Us_s) Cf_Ma(Us_s,U,vs), Us0);
    
    function z0 = Cf_Ma(Us_s, U, vs)
        sus = Us_s./vs;
        z0 = log10(Us_s.^2/U.^2)+1.53*(log10(sus)).^2+0.15*log10(sus)+1.91;
    end

end