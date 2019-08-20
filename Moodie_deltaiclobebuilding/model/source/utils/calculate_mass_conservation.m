function [mc, calc_idx] = calculate_mass_conservation(s, Bc, Be, Bc0, Bf0, Bo0, rad_idx, sig, lastavul_idx, dx, con)
    %% calculate the volume of lobe and delta topset regions
    % this establishes a mass conservation
    
    %%%
    %%% NOTE: these routines only works up to a single avulsion cycle!!
    %%%
    
    nt = size(s.eta, 2);
    skip = 500;
    cnt = 1;
    calc_idx = [1:skip:nt, nt];
    nx = size(s.eta, 1);
    
    
    
    
    for t = calc_idx
        disp(['t: ', num2str(t), ' of ' num2str(nt)])
        
        %% discretized calculations
        deta_lastavul = (  s.eta(:, t)'  ) - (  s.eta(:, lastavul_idx) - (sig*t*86400) )';  

        % volume in delta topset region including in channel
        vol_in_rad_disc = ( deta_lastavul(1:rad_idx) .* [dx/2 repmat(dx, 1, rad_idx-2) dx/2] .* Be(1:rad_idx) ) .* (1-con.phi);
        vol_in_rad_disc_t = sum(vol_in_rad_disc);

        % volume in the floodplain for redist
        vol_floodplain_disc = ( deta_lastavul(1:rad_idx) .* [dx/2 repmat(dx, 1, rad_idx-2) dx/2] .* (Be(1:rad_idx)-Bc(1:rad_idx)) ) .* (1-con.phi);
        vol_floodplain_disc_t = sum(vol_floodplain_disc);
        
        % volume in the lobe for redist (method 1)
        detadt = horzcat( zeros(size(s.eta(:,1))), (s.eta(:, 2:t) - s.eta(:, 1:t-1)) ); % change in elecation at every space through time to t
        timevar_lobe = sum(detadt(rad_idx:end, :) .* s.Be(rad_idx:end, 1:t), 2); % getadt times the depositional width at time of deta
        timevar_Vol_lobe_x = timevar_lobe' .* [dx/2 repmat(dx, 1, length(detadt(rad_idx:end, 1))-2) dx/2] * (1-con.phi); % multiplied by the constant space vector
        vol_to_lobe_disc_t = sum(timevar_Vol_lobe_x); % cumulate
        
        sum ( sum( detadt(rad_idx:end, :) .* s.Be(rad_idx:end, 1:t) .* repmat([dx/2 repmat(dx, 1, length(detadt(rad_idx:end, 1))-2) dx/2]', 1, t), 2) );
        sum(  trapz(detadt(rad_idx:end, :) .* s.Be(rad_idx:end, 1:t),2)*dx*(1-con.phi)  );
        
        % volume in the lobe for redist (method 2)
        vol_to_lobe_disc2 = ( deta_lastavul(rad_idx:end) .* [dx/2 repmat(dx, 1, nx-rad_idx)] .* Be(rad_idx:end) ) .* (1-con.phi);
        vol_to_lobe_disc_t2 = sum(vol_to_lobe_disc2);
        
        
        %% vector tracked calculations
        
        % volume input to the model domain
        idx = 1;
        vol_input_vect = integrate_qs_flux(s.qs(:, 1:t), s.Bc(:, 1:t), idx);
        
        % volume in lobe region
        idx = rad_idx;
        vol_to_lobe_vect_t = integrate_qs_flux(s.qs(:, 1:t), s.Bc(:, 1:t), idx);
        
        % volume stored in radius
        vol_in_rad_vect_t =  vol_input_vect - vol_to_lobe_vect_t;
        
        % volume stored in floodplain
        vol_floodplain_vect_t = vol_in_rad_vect_t .* (Bf0 / (Bc0+Bf0));
    
        
        %% store values
        mc.vol_to_lobe_disc(cnt) = vol_to_lobe_disc_t;
        mc.vol_to_lobe_vect(cnt) = vol_to_lobe_vect_t;
        mc.vol_in_rad_disc(cnt) = vol_in_rad_disc_t; 
        mc.vol_in_rad_vect(cnt) = vol_in_rad_vect_t;
        mc.vol_floodplain_disc(cnt) = vol_floodplain_disc_t;
        mc.vol_floodplain_vect(cnt) = vol_floodplain_vect_t;
        mc.vol_input_vect(cnt) = vol_input_vect;
        
        cnt = cnt + 1;
        
        clear vol_to_lobe_disc_t vol_to_lobe_vect_t vol_in_rad_disc_t vol_in_rad_vect_t vol_floodplain_disc_t vol_floodplain_vect_t vol_input_vect

    end
    
end