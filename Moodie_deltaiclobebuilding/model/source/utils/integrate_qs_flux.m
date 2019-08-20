function [volume_int] = integrate_qs_flux(qsmat, Bcmat, idx)
    % takes two matricies subsetted from s structure to calcualte the flux
    % of sediment across the plane defined at idx
    
    volume_int = sum( trapz(qsmat(idx, :) .* double(Bcmat(idx, :)))  * 86400);

end