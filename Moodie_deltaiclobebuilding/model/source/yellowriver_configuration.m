function [cfg] = yellowriver_configuration()
    %% configuration parameters that should be constant for all model runs relating to the Yellow River

	% constants
    cfg.con.g = 9.81;                           % gravitational constant
    cfg.con.rho_f = 1000;                       % fluid density
    cfg.con.rho_s = 2650;                       % sediment density
    cfg.con.Rr = (cfg.con.rho_s / ...
                        cfg.con.rho_f) - 1;     % submerged specific gravity
    cfg.con.phi = 0.4;                          % porosity
    cfg.con.nu = 1.004 * 1e-6;                  % kinematic viscosity, m^2/s
    
    % spatial domain params
    cfg.L = 400e3;                              % total model length, meters
    cfg.nx = 600;                               % number of spatial steps, [600]
    cfg.shoreleni = 80e3;                       % initial length of shoreline, determines delta apex
    
    % temporal domain params
    cfg.dti = 0.25;                             % initial timestep, explicit to fraction of day
    cfg.dtsec = cfg.dti * 86400;                % timestep initial in seconds
    cfg.dtsec_int = 86400 .* [0.0001, 0.5];     % limits to the timestep: [0.00005, 0.5] is a stable choice
    
    % stability params
    cfg.au = 0.8;                               % winding coefficient
    cfg.CFLc = 0.00005;                         % timestep computation adjustment
    cfg.cfl_idx_fun = (@(mou_idx) mou_idx - 1); % index to calculate cfl at
    cfg.psi = 25;                               % foreset slope maximum
    
    % real-world params
    cfg.D50 = 90e-6;                            % D50 in meters
    cfg.S0 = 6.4e-5;                            % initial fan slope
    cfg.sigma = 5;                              % subsidence rate, mm/yr
    cfg.Cf = 0.001;                             % friction coefficient (backwater and sediment transport)
    cfg.Qwbf = 3000;                            % bankfull discharge, m3/s
    cfg.Qwform = 1300;                          % formative discharge, m3/s (progradation crit.)
    cfg.H0 = 0;                                 % downstream boundary for backwater calc
    
    % sediment transport params
    [p1, p2] = set_GEH_params('HM');            % coefficients to Generalized Engelund Hansen formulation
    cfg.p1 = p1;
    cfg.p2 = p2;

    % width-related params
    cfg.Bc0 = 400;                              % channel width
    cfg.Bf0 = cfg.Bc0*10;                       % floodplain width
    cfg.Bo0 = 9e3;                              % overbank width
    cfg.thet = 5;                               % plume spreading angle
    cfg.gamm = 90;                              % fan delta angle, degrees
    
    % other parameters
    cfg.preavul_close_thresh = 0.00;            % values within this percentage are allowed to avulse instead

end