function [s, analysis] = virtualdelta(cfg)


    %% unpack the cfg structure here
    % see configuration files for variable definitions
    con = cfg.con;
    S0 = cfg.S0;
    L = cfg.L;
    nx = cfg.nx;
    cfg.dx = L / cfg.nx;
    dx = cfg.dx;
    x = 0 : dx : L;
    
    au = cfg.au;
    CFLc = cfg.CFLc;
    dtsec = cfg.dtsec;
    dtsec_int = cfg.dtsec_int;
    cfl_idx_fun = cfg.cfl_idx_fun;
    psi = cfg.psi;
    
    Bc0 = cfg.Bc0;
    Bf0 = cfg.Bf0;
    Bo0 = cfg.Bo0;
    Qwbf = cfg.Qwbf;
    Qwform = cfg.Qwform;
    thet = cfg.thet;
    gamm = cfg.gamm;
    
    Cf = cfg.Cf;
    H0 = cfg.H0;
    
    D50 = cfg.D50;
    p1 = cfg.p1;
    p2 = cfg.p2;
    sig = cfg.sigma / 1000 / 31536e3; % subsidence coverted to m/s
    Int = cfg.Int;
    
    storage = cfg.storage;
    
    mouth_switch = cfg.mouth_switch;
    preavul_switch = cfg.preavul_switch;
    preavul_thresh = cfg.preavul_thresh;
    preavul_trigg = cfg.preavul_trigg;
    preavul_close_thresh = cfg.preavul_close_thresh;
   
    % convert nan input to ignore overbank flow
    if and(isnan(preavul_trigg), ~(strcmp(preavul_switch, 'qingshuigou')))  
        preavul_switch = 'any-setup+topset';
        warning('preavulSwitch changed due to NaN value')
    end
    
    %% model domain setup
    Hnbf = Hn_DW(Qwbf, Bc0, Cf, S0, con.g); % bankfull flow depth (i.e., formative discharge flow depth)
    Hnform = Hn_DW(Qwform, Bc0, Cf, S0, con.g); % formative depth
    Lblow = Hnform / S0;
    
    Qw(:, 1) = set_Qw(cfg.Qw_str, cfg.Qw_num);
    
    % use configuration param to setup initial channel
    zed_zero_idx = floor(nx/2);
    [eta, zed, rad] = setup_initial_channel(cfg.initial_channel, zed_zero_idx, Hnbf, S0, ...
                                            nx, dx, repmat(Bc0, 1, nx+1), H0, Cf, Qwbf);

    eta_yi = eta(1); % used for avulsion calculations
    
    radi = rad;
    moui = radi; % inital mouth location, meters
    mou_idx = get_idx(x, moui);
    rad_idx = get_idx(x, radi);
    mou0_idx = mou_idx;
    
    gammleni = cfg.shoreleni / (gamm/180) / pi; % initial radius of the delta from rad to apex
    gammapex = rad - gammleni;
    gammapex_idx = get_idx(x, gammapex);
    shore_leni = (1/4) * (2 * pi * (radi - gammapex)) / 1000; % deg2rad(gamm)
    [Bc] = set_Bc(mou_idx, Bc0, thet, nx, dx); % Bc is channel width, Be is Exner width
    [Be] = set_Be(Bc0, Bf0, Bo0, mou_idx, rad_idx, nx); % Bc is channel width, Be is Exner width
    
    % time looping counters
    cnt.l = 1; % time loop counter
    cnt.P = 2; % print counter
    cnt.y = 0; % year counter
    cnt.d = 1; % day counter
    cnt.dfrac = 0; % day fraction counter
    cnt.delap = 0; % count days elapsed, no resets
    cnt.nAvuls = 0; % number of avulsions, no resets
    cnt.dsecs = 0; % day seconds elapsed
    cnt.Qw_num = cfg.Qw_num; % counter for moving Qw track along

    % set up the run_ctrl object
    run_ctrl = run_ctrl_class(cfg.run_str, cfg.run_n);
    
    % set some conditions
    eta_i = eta; % save initial eta profile
    s.eta(:, 1) = eta_i; % save initial eta profile
    s.mou_idx(1) = mou_idx;
    s.rad(1) = rad;
    flux_to_lobe = 0;
    flux_in_rad = 0;
    lastavul_idx = 1; % one before initial time
    avulloc = 1;
    eta_lastavul = eta;
    avulsion_this_year = false;
    
    % add some more elements to storage for drawing
    d.x = x;
    d.dx = dx;
    d.Lblow = Lblow;
    d.Hnform = Hnform;
    d.Hnbf = Hnbf;
    d.Bf0 = Bf0;
    d.Bc0 = Bc0;
    d.Bo0 = Bo0;
    d.gamm = gamm;
    d.gammapex = gammapex;
    d.lij_loc = moui - 56e3; % location of lijin
    [d.lij_idx] = get_idx(x, moui - 56e3);
    
    if cfg.runtime_plot
        % setup runtime figure if configured
        runtime_fig = figure(); 
        hold on;
        cnt.runtime_plot = 0;
    end
    
    
    while run_ctrl.status()
        
        % calculate bed slope
        [S] = calculate_bed_slope(eta, nx, dx);
        
        % calculate backwater
        [H] = calculate_backwater_dBdx(eta, S, Bc, H0, Cf, Qw(cnt.delap+1), nx, dx);
        
        % backwater location
        [Lb] = find(((H(2:end) - H(1:end-1)) / dx) > 5e-6, 1, 'first');
        
        % caclulate velocity
        [U] = calculate_velocity(Qw(cnt.delap+1), H, Bc);
        
        % calculate sediment transport
        [qs] = calculate_qs_GEH(U, Cf, D50, p1, p2, con);
        qsBC = qs(1)*1.0005; % upstream boundary condition
        
        % calculate new timestep based on sediment transport
        if 86400 - cnt.dsecs < dtsec
            dtsec = 86400 - cnt.dsecs;
        else
            cfl_idx = cfl_idx_fun(mou_idx)+4;
            [dtsec] = calculate_new_dtsec(qs(cfl_idx)-qs(cfl_idx-1), Be(cfl_idx), Bc(cfl_idx), dx, CFLc, con, dtsec_int);
        end
        
        % record flux of sediment to the lobe
        flux_to_lobe_dt = ( qs(rad_idx) * dtsec * Bc0);
        flux_to_lobe = flux_to_lobe + flux_to_lobe_dt;
        flux_in_rad =  flux_in_rad + (( qs(1) * dtsec * Bc0) - flux_to_lobe_dt);
        
        % calculate new bed elevation
        [eta] = calculate_new_bed_exner(eta, au, qs, qsBC, Bc, Be, nx, dx, dtsec, Int, con);
        
        % apply subsidence
        [eta] = apply_subsidence(eta, zed, sig, dtsec);
        
        % calculate deposit thickness/attributes
        [dep] = get_deposit_thick(eta, eta_lastavul);
        
        % avalanche foreset if needed
        [eta] = check_for_avalanche(eta, S, Be, mou_idx, dep, psi);
        
        % check for progradation of lobe
        [mou_idx] = check_for_mou_prog(eta, zed, dep, Hnform, rad_idx, mou_idx, mouth_switch);
        if mou_idx ~= mou0_idx
            % if mouth location has changed, update the model configs
            [Bc] = set_Bc(mou_idx, Bc0, thet, nx, dx);
            [Be] = set_Be(Bc0, Bf0, Bo0, mou_idx, rad_idx, nx);
            mou0_idx = mou_idx;
        end
        
        % add a time/loop count
        cnt.l = cnt.l + 1;
        cnt.dsecs = cnt.dsecs + dtsec;
        
        % if equal to one day has passed
        if cnt.dsecs >= 86400
           
            if cfg.runtime_plot && mod(cnt.delap, cfg.runtime_plot_int) == 0
                % temporary debugging plot
                [runtime_fig] = runtime_plot(runtime_fig, s, x, eta, zed, H, Hnform, U, Qw, mou_idx, rad_idx, Bc, Be, dx, cnt);

                if cfg.runtime_save
                    pause(0.1)
                    set(gcf, 'Pos', [50 100 900 500], 'PaperPositionMode', 'auto')
                    print('-dpng', '-r300', sprintf('./movies/runtime/%04d.png', cnt.runtime_plot))
                end
                
                cnt.runtime_plot = cnt.runtime_plot + 1;
                
            end
            
            % stdout updater
            if mod(cnt.delap, cfg.runtime_stdout_int) == 0
                cnt.n = run_ctrl.counter;
                cnt.N = run_ctrl.limit;
                cnt.P = cnt.P + 1;
                stdout_updater(cnt)
            end
            
            % handle counting for the daily varying Qw
            cnt.d = cnt.d + 1; % add a daystep
            cnt.delap = cnt.delap + 1; % add an elapsed day
            cnt.dfrac = dtsec / 86400; % reset day fraction counter
            cnt.dsecs = 0;
            
            if run_ctrl.day_incr
                run_ctrl = run_ctrl.increment(); % increment the run_ctrl counter
            end
            
            % check avulsion criteria in the model
            [activate, avulloc] = check_for_avulsion(eta, zed, H, Hnbf, rad_idx, gammapex_idx, ...
                preavul_switch, preavul_thresh, preavul_trigg, preavul_close_thresh, nx);
            
            % if avulsion criteria is met and no avulsion this flood year
            if and(activate, ~avulsion_this_year)
                
                % calculate time since the last avulsion (seconds)
                time_since = (cnt.delap - lastavul_idx) * 86400;
                
                % calculate volume of sed in lobe
                [Vol_floodplain, Vol_lobe] = calculate_avulsion_volumes(eta, Bc0, Bf0, Bc, Be, rad_idx, eta_lastavul, flux_in_rad, flux_to_lobe, time_since, sig, dx, con);
                
                % complete avulsion in planform domain
                rad0 = rad;
                zed0 = zed;
                [rad, zed] = avulsion_plan(rad0, zed0, gamm, gammapex, gammapex_idx, Bf0, S0, eta_yi, Vol_floodplain, Vol_lobe, x, dx, con);
                [rad_idx] = get_idx(x, rad);
                mou_idx = rad_idx;
                
                % complete avulsion in long profile
                eta0 = eta;
                [eta] = avulsion_long(eta_i, eta0, zed, rad_idx, Hnbf, avulloc, dx, cfg.postavul_switch);
                [Bc] = set_Bc(mou_idx, Bc0, thet, nx, dx);
                [Be] = set_Be(Bc0, Bf0, Bo0, mou_idx, rad_idx, nx);
                
                % reset/update params
                lastavul_idx = cnt.delap;
                cnt.nAvuls = cnt.nAvuls + 1;
                eta_lastavul = eta;
                flux_to_lobe = 0;
                flux_in_rad = 0;
                avulsion_this_year = true;
                
                if run_ctrl.avul_incr
                    run_ctrl = run_ctrl.increment(); % increment the run_ctrl counter
                end
            else
                % if not activated, need to reset avulloc to zero so nothing is logged
                avulloc = NaN;
            end
            
            % store data for end of run
            if storage
                [s] = data_to_storage(s, cnt, eta, H, zed, mou_idx, avulloc, rad, Qw(cnt.delap), Lb, Be, Bc, qs);
            else
                [s] = data_to_storage_lite(s, cnt, mou_idx, avulloc, rad);
            end
            
            % check if a year has passed to reset Qw curves and avulsion_this_year
            if cnt.d == 366 % if a year has passed
                cnt.y = cnt.y + 1; % add a yearstep
                cnt.d = 1; % reset day counter
                cnt.Qw_num = run_ctrl.Qw_func(cnt.Qw_num);
                Qw(:, end+1) = set_Qw(cfg.Qw_str, cnt.Qw_num);
                avulsion_this_year = false;
            end
            
        % end else (day checker) statement    
        end
        
        if any(eta > 100) % || any(isnan(eta))
            run_ctrl = run_ctrl.kill();
            disp('Error, unstable: produced eta > 100')
            disp(['cfl_idx', num2str(cfl_idx)])
            disp(['unstable idx', num2str(find(eta>100))])
        end
        
    % end time while loop    
    end
    
    stdout_updater(cnt)
    
    % save the raw data
    if ~isnan(cfg.savedata_filename) % if savedata_filename name is given
        all_vars = whos;
        vars_save = cellfun(@isempty, regexp({all_vars.class}, '^matlab\.(ui|graphics)\.')); % variable names to save (no figures)
        save(cfg.savedata_filename, all_vars(vars_save).name, '-v7.3')
    end
    
    % save the final just-before-avulsion-state (for engineered avulsion spinups)
    if ~isnan(cfg.savefinal_filename) % if savedata_filename name is given
        all_vars = whos;
        vars_save = cellfun(@isempty, regexp({all_vars.class}, '^matlab\.(ui|graphics)\.')); % variable names to save (no figures)
        no_s = ~strcmp({all_vars.name}, 's');
        save(cfg.savefinal_filename, all_vars(and(vars_save, no_s)).name, '-v7.3')
    end
    
    % do computations/plotting
    if ~isnan(cfg.analysis_filename) % if analysis_filename is given
        [d] = extract_avulsion_information(s, d);
        [d] = extract_additional_information(s, d);
        [analysis] = convert_to_analysis(s, d);
        analysis.cfg = cfg;
        save(cfg.analysis_filename, 'analysis');
    end
    
    s.d = d;
    
    % s gets returned as 'ans', in function definition

% end function definition for virtualdelta    
end



