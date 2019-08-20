function [eta] = avulsion_long(etai, eta0, zed, rad_idx, Hnbf, avulloc, dx, postavulSwitch)

    win = 21; % width of taper, i.e. number of nodes from loc to taper
    hwin = floor(win/2);
    
    Hnbf_cut = zed - Hnbf;
    
    if isnan(postavulSwitch)
        error('NaN given as postavulSwitch, no avulsions permitted, bad config')
    end
    
    switch postavulSwitch
        case 'initial'
            eta = etai;
        case 'taper-to-etai'
            % taper down to the initial bed topography beyond avulsion
            eta = eta0; % make new bed vector
            [taper, taper_idx] = make_bed_taper(eta0, etai, avulloc-hwin, win, dx); % calculate taper
            eta(taper_idx) = taper;
            taper_toe_idx = find(taper_idx, 1, 'last');
            eta(taper_toe_idx:end) = etai(taper_toe_idx:end);
        case 'taper-to-bankfull'
            % taper down to a bakfull flow depth below zed out to rad
            eta = eta0; % make new bed vector
            datum = Hnbf_cut;
            [taper, taper_idx] = make_bed_taper(eta0, datum, avulloc-hwin, win, dx); % calculate taper
            eta(taper_idx) = taper;
            taper_toe_idx = find(taper_idx, 1, 'last');
            
            [runout, runout_idx] = make_bed_runout(eta, etai, Hnbf_cut, taper_toe_idx, rad_idx); 
            eta(runout_idx) = runout;

            mesh_idx = taper_toe_idx:find(eta < taper(end), 1, 'first');
            eta(mesh_idx) = linspace(taper(end), eta(mesh_idx(end)), length(mesh_idx));
            
    end
end