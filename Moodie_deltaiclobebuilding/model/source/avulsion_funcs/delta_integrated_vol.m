function [vol] = delta_integrated_vol(rad, eta0, S0, gamm, gammapex)
    %% volumetric calculation for "delta volume" for given rad
    % solve for the volume for a given `rad` based on the configuration of
    % gammapex. This uses volume of cylinder and cone formulas to
    % difference out the volume of the region of interest by revolving
    % about the gammapex
    %
    frac_circ = (deg2rad(gamm)/(2*pi)); % fraction of circle the delta angle represents
    x_zero = eta0 / S0; % x coordinate of y=0
    h_lower = abs(rad*-S0+eta0);
    h_upper = (gammapex*-S0+eta0) - (x_zero*-S0+eta0); % eta at gammapex - eta at x_zero
    upper_cone = (pi*frac_circ) .* (x_zero-gammapex).^2 .* ( h_upper / 3 ); % upper cone between x_zero and gammapex
    main_cylinder = (pi*frac_circ) * (rad-gammapex).^2 .* ( h_lower ); % main cylinder producing the positive delta volume
    lower_cone = (pi*frac_circ) .* (rad-gammapex).^2 .* ( h_lower / 3 ); % lower cone representing the basin bottom, or antecedent topo
    vol = main_cylinder - (lower_cone - upper_cone);
end