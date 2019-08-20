function [rad, zed] = avulsion_plan(rad0, zed0, gamm, gammapex, gammapex_idx, Bf0, S0, eta0, Vol_floodplain, Vol_lobe, x, dx, con)
    %% execute the avulsion routine "planform" parts
    % this includes the progradation of the delta front and the aggradation
    % of the delta topset.
    
    %% calculate and prograde based on the lobe volume
    rad0_idx = get_idx(x, rad0); % old radius index
    shorelen = (rad0 - gammapex) * deg2rad(gamm); % length of the shoreline at old radius
    
    %% some hacky numerical integration of the volume of delta increase for
    % every 1 m increase in radius. This is used to determine how much the
    % volume of the lobe should cause the delta radius to increase by.
    idx_lim = 20000; % limit (m) to calculate to from the current radius; 
    
    drad = [];
    while isempty(drad)
        % integration, uses distance, basin depth (slope), and shore length
        %
        [vol0] = delta_integrated_vol(rad0, eta0, S0, gamm, gammapex); % current radius delta volume

        % define coordinates to evaluate new volume at, and evaluate
        int_coords = 0:1:idx_lim;
        for i = 1:length(int_coords)
            [voln(i)] = delta_integrated_vol(rad0+int_coords(i), eta0, S0, gamm, gammapex);
        end

        % calculate change in volume with every meter
        dvol = [voln(1)-vol0, voln(2:end) - voln(1:end-1)];
        cumdvol = cumsum(dvol);

        % find closest match
        drad = find( (sum(Vol_lobe)./(1-con.phi)) < cumdvol, 1, 'first' ); % find the best match
        % do an interpolation?
        
        idx_lim = floor(idx_lim * 1.2);
    end
    
    % new radius is the old plus change
    rad = rad0 + drad; % new radius
    rad_idx = get_idx(x, rad);

    %% calculate and aggrade based on the topset volume
    % set up a master vector for the area corresponding to any given x
    agg_coords = [x(2:end)-dx/2, x(end)+dx/2]; % coordinates offset from x to do area calculations at
    agg_idx = 1:rad0_idx; % index of coordinates which should be eligable to be aggraded
    A_agg_coords = ones(size(x)); % preallocate a vector to be filled with depositional areas
    
    % width in the neck is the constant floodplain width
    A_agg_coords(1:gammapex_idx+1) = Bf0*dx; % width in the neck
    
    % width in delta is like a pie slice
    pie_coords = agg_coords(gammapex_idx+1:rad0_idx+1) - gammapex;
    pie_area_fun = (@(rad_, gamm) (pi .* rad_.^2) .* (deg2rad(gamm)/(2*pi)));
    pie_area_cum = pie_area_fun(pie_coords, gamm);
    pie_area_slice = [0, pie_area_cum(2:end) - pie_area_cum(1:end-1)];
    pie_area_slice(pie_area_slice < Bf0*dx) = Bf0*dx; % nothing is less than neck width
    A_agg_coords(gammapex_idx:rad0_idx) = pie_area_slice; % replace delta region with pie areas
    Vol_floodplain_padded = [Vol_floodplain, zeros(1, length(x)-length(Vol_floodplain))];
    zed = zed0 + (   Vol_floodplain_padded ./ (1-con.phi) ./ (A_agg_coords)   ); % elementwise aggradation
    
end
