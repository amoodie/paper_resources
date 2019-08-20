function [eta, zed, rad] = setup_initial_channel(initial_channel, zed_zero_idx, Hnbf, S0, nx, dx, Bc, H0, Cf, Qwbf)

    if strcmp(initial_channel, 'new')
        % setup a new channel based on the config given in the cfg
        [zed] = set_zedi(zed_zero_idx, S0, nx, dx);
        [eta] = set_etai(zed, Hnbf, S0, dx);

        % find radius
        rad = find(zed <= 0, 1)*dx; % initial delta radius (from model boundary), meters, intersection of delta with sea level
        
        % make adjustments to zed, add params
        zed(zed<0) = 0;
    else
        error('NotImplementedError')
        
    end
end