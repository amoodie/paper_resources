function [caseStruct] = MellorYamadaConcVelModel(D_dim, cb_i, u_star, S_0, h)
%% WPPRedictConcVelModel makes predictions for concentration and velocity profiles
%
%   This is the main function, which should be given the stationStruct object,
%   and returns a new stationStruct with the attributes of stratified and unstratified models
%   
%   The methods rely HEAVILY on the work by Yeh and Parker (https://doi.org/10.1016/j.cageo.2011.12.004)
%
   
    % sanity checks
    if ~( length(D_dim) == length(cb_i) )
        error('input lengths must match (number of classes)')
    end
    
    if any([ length(u_star)>1, length(S_0)>1, length(h)>1 ])
        error('u_*, S_0, and h should be scalars')
    end
    
    % rotate so that all vectors are in correct orientation (row vectors)
    if ~isrow(D_dim)
        D_dim = D_dim';
    end
    
    if ~isrow(cb_i)
        cb_i = cb_i';
    end
    
    % a few more params
    D_frac = length(D_dim); % number of grain size fractions
    C_frac = cb_i / sum(cb_i); % the fractions in each class
    C_bed_frac = C_frac; % copy to bed, just for clarity of method in a few points of code
    
    b = 0.05;
    Nz = 51;
    
    %% set up grain size data for conveniece here
    [gsMetrics] = calculate_grainsize_metrics(u_star, D_dim, C_bed_frac);


    %% do the simple (log and Rouse) models
    [simple] = simpleModel(u_star, C_bed_frac, gsMetrics.V, S_0, h, b, Nz);
    
    
    %% do the full model with concentration == 0
    solution_plot = false; 
    final_plot = false;
    verbose = false;
    [strat_nosed] = stratModel(u_star, 0, 0, 1, 0, S_0, b, h, Nz, solution_plot, final_plot, verbose, simple);
    
    
    %% do the full model with real values
    solution_plot = false; 
    final_plot = false;
    verbose = false;
    [strat] = stratModel(u_star, D_dim, cb_i, C_bed_frac, gsMetrics.V, S_0, b, h, Nz, solution_plot, final_plot, verbose, simple);

    
    %% do some computations and stick data into the stationStruct
    caseStruct.nu_t_non       = simple.nu_t_non .* (u_star*h);
    caseStruct.nu_t           = strat.nu_t .* (u_star*h);
    caseStruct.nu_t_nosed     = strat_nosed.nu_t .* (u_star*h);

    caseStruct.Us_non         = simple.u_non      .* u_star;
    caseStruct.Us             = strat.u_hat       .* u_star;
    caseStruct.Us_nosed       = strat_nosed.u_hat .* u_star;

    caseStruct.CsSum_non     = simple.c_non      .* sum(cb_i);
    caseStruct.CsSum         = strat.CsSum;
    caseStruct.CsSum_nosed   = strat_nosed.CsSum;

    caseStruct.Cs_non        = simple.c_i_non      .* sum(cb_i);
    caseStruct.Cs            = strat.Cs;
    caseStruct.Cs_nosed      = strat_nosed.Cs;

    caseStruct.k_hat_non     = NaN(size(strat.k_hat));
    caseStruct.k_hat         = strat.k_hat;
    caseStruct.k_hat_nosed   = strat_nosed.k_hat;

    caseStruct.k_non         = NaN(size(strat.k));
    caseStruct.k             = strat.k;
    caseStruct.k_nosed       = strat_nosed.k;

    caseStruct.l_0_non       = NaN(size(strat.k));
    caseStruct.l_0           = strat.l_0;
    caseStruct.l_0_nosed     = strat_nosed.l_0;
    
    caseStruct.Zs_non       = simple.z;
    caseStruct.Zs           = strat.Zs;
    caseStruct.Zs_nosed     = strat_nosed.Zs;
    
    caseStruct.K_red = mean(caseStruct.nu_t) / mean(caseStruct.nu_t_nosed);

    
    %% make a plot
    combine_plot = false;
    if combine_plot
        close all
        figure()
        subplot(1,4,1); hold on;
            cla
            plot(simple.u_non .* u_star, simple.eta, 'k--')
            plot(strat.u_hat .* u_star, strat.eta, 'r')
            plot(strat_nosed.u_hat .* u_star, strat_nosed.eta, 'b')
            plot(caseStruct.Velocity.data.velMean, caseStruct.Velocity.data.measZ./h, 'ko')
            xlabel('$u$ (m/s)')
            ylabel('$z/h$')
        subplot(1,4,2); hold on;
            cla
            plot(simple.c_non.*sum(cb_i), simple.eta, 'k--')
            plot(strat.c_hat.*strat.C_mean, strat.eta, 'r')
            plot(strat_nosed.c_hat.*strat_nosed.C_mean, strat_nosed.eta, 'b')
            plot(caseStruct.waterSamplesTable.concNW./2650, caseStruct.waterSamplesTable.sampleZ./h, 'ko')
            xlabel('$c$ (-)')
        subplot(1,4,3); hold on
            cla
            plot(simple.nu_t_non .* (u_star*h), simple.eta, 'k--')
            plot(strat.nu_t .* (u_star*h), strat.eta, 'r')
            plot(strat_nosed.nu_t .* (u_star*h), strat_nosed.eta, 'b')
            xlabel('$\nu_t$ (m$^2$/s)')
        subplot(1,4,4); hold on;
            cla
            plot(strat.k_hat, strat.eta, 'r')
            plot(strat_nosed.k_hat, strat_nosed.eta, 'b')
            xlabel('$k/u_*^2$ (m$^2$/s$^2$)')
        drawnow 
    end
   
    
        
end
