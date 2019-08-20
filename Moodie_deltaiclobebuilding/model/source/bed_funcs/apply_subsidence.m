function [eta, zed] = apply_subsidence(eta, zed, sig, dtsec)
    % suside eta and floodplain
    % sig is passed in units of m/s
    
    eta0 = eta;
    zed0 = zed;
    dz = -sig * dtsec;
    eta = eta0 + dz;
    zed = zed0 + dz;
end