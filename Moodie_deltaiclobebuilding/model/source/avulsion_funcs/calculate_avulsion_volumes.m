function [Vol_floodplain, Vol_lobe, Vol_inside_rad] = calculate_avulsion_volumes(eta, Bc0, Bf0, Bc, Be, rad_idx, eta_lastavul, flux_in_rad, flux_to_lobe, time_since, sig, dx, con)
    %% calculate the volume of lobe and delta topset regions
    % this is needed for redistribution of volume following avulsion
    % time since  == time since last avulsion in seconds
    
    % change in eta since last avulsion
    deta_lastavul = eta - (eta_lastavul - (sig * time_since)); % shift the last surface down to account for subsidence
    
    % volume in delta topset region including in channel
    Vol_inside_rad = ( deta_lastavul(1:rad_idx) .* [dx/2 repmat(dx, 1, rad_idx-2) dx/2] .* Be(1:rad_idx) ) .* (1-con.phi);
    
    % volume in the floodplain for redist
    Vol_floodplain_disc = ( deta_lastavul(1:rad_idx) .* [dx/2 repmat(dx, 1, rad_idx-2) dx/2] .* (Be(1:rad_idx)-Bc(1:rad_idx)) ) .* (1-con.phi);
    
    % determine error in mass balance calculation for x-coordinate specific calculation
    % error is due to discretized spatial calculation
    Vol_floodplain_flux = flux_in_rad * (Bf0 / (Bc0+Bf0));
    mass_balance_perc_error = (sum(Vol_floodplain_disc) - Vol_floodplain_flux) / Vol_floodplain_flux * 100;
    mass_balance_error = sum(Vol_floodplain_disc) - Vol_floodplain_flux;

    % correct the x-coord specific one for the mass error, spreading evenely
    Vol_floodplain = Vol_floodplain_disc + (-1 .* mass_balance_error ./ length(Vol_floodplain_disc));
    
    % volume in lobe region
    Vol_lobe = flux_to_lobe;
    
end