function [stationStruct] = MellorYamadaConcVelModel_struct(stationStruct)
%% WPPRedictConcVelModel makes predictions for concentration and velocity profiles
%
%   This is the main function, which should be given the stationStruct object,
%   and returns a new stationStruct with the attributes of stratified and unstratified models
%   
%   The methods rely HEAVILY on the work by Yeh and Parker (https://doi.org/10.1016/j.cageo.2011.12.004)
%
   
    %% extract some real variables for evaluating the models
    cb_i = stationStruct.concProf.gsClassPred.params.cb / 2650;
    u_star = stationStruct.Velocity.ustarCalib; % old format, load it form stationStruct (g*h*S)^(0.5);

    S_0 = stationStruct.Velocity.slope; 
    h = stationStruct.FlowDepthAtCollection;

    gsTable = stationStruct.nearBedData.gsDistNearBedNWnorm;
    D_frac = size(gsTable, 1); % number of grain size fractions
    C_frac = gsTable.avgClassNorm' ./ 100; % the fractions in each class
    C_bed_frac = C_frac; % copy to bed, just for clarity of method in a few points of code
    D_dim = (gsDistClass2num(gsTable) .* 1e-6)'; % dimensions grain sizes, (m)
    
    b = 0.05;
    Nz = 51;
    
    %% set up grain size data for conveniece here
    [gsMetrics] = calculate_grainsize_metrics(u_star, D_dim, C_bed_frac);


    %% do the simple (log and Rouse) models
    [simple] = simpleModel(u_star, C_bed_frac, gsMetrics.V, S_0, h, b, Nz);
    
    
    %% do the full model with concentration == 0
    solution_plot = false; 
    final_plot = false;
    verbose = true;
    disp( ['            simulating without sediment....'] )
    [strat_nosed] = stratModel(u_star, 0, 0, 1, 0, S_0, b, h, Nz, solution_plot, final_plot, verbose, simple);
    
    %% do the full model with real values
    solution_plot = false; 
    final_plot = false;
    verbose = true;
    disp( ['            simulating with sediment....'] )
    [strat] = stratModel(u_star, D_dim, cb_i, C_bed_frac, gsMetrics.V, S_0, b, h, Nz, solution_plot, final_plot, verbose, simple);

    
    %% do some computations and stick data into the stationStruct
    stationStruct.velProf.MYfullPred.nu_t_non       = simple.nu_t_non .* (u_star*h);
    stationStruct.velProf.MYfullPred.nu_t           = strat.nu_t .* (u_star*h);
    stationStruct.velProf.MYfullPred.nu_t_nosed     = strat_nosed.nu_t .* (u_star*h);

    stationStruct.velProf.MYfullPred.Us_non         = simple.u_non      .* u_star;
    stationStruct.velProf.MYfullPred.Us             = strat.u_hat       .* u_star;
    stationStruct.velProf.MYfullPred.Us_nosed       = strat_nosed.u_hat .* u_star;

    stationStruct.concProf.MYfullPred.CsSum_non     = simple.c_non      .* sum(cb_i)          .* 2650;
    stationStruct.concProf.MYfullPred.CsSum         = strat.CsSum;
    stationStruct.concProf.MYfullPred.CsSum_nosed   = strat_nosed.CsSum;

    stationStruct.concProf.MYfullPred.Cs_non        = simple.c_i_non      .* sum(cb_i)          .* 2650;
    stationStruct.concProf.MYfullPred.Cs            = strat.Cs;
    stationStruct.concProf.MYfullPred.Cs_nosed      = strat_nosed.Cs;

    stationStruct.concProf.MYfullPred.k_hat_non     = NaN(size(strat.k_hat));
    stationStruct.concProf.MYfullPred.k_hat         = strat.k_hat;
    stationStruct.concProf.MYfullPred.k_hat_nosed   = strat_nosed.k_hat;

    stationStruct.concProf.MYfullPred.k_non         = NaN(size(strat.k));
    stationStruct.concProf.MYfullPred.k             = strat.k;
    stationStruct.concProf.MYfullPred.k_nosed       = strat_nosed.k;

    stationStruct.concProf.MYfullPred.l_0_non       = NaN(size(strat.k));
    stationStruct.concProf.MYfullPred.l_0           = strat.l_0;
    stationStruct.concProf.MYfullPred.l_0_nosed     = strat_nosed.l_0;
    
    stationStruct.concProf.MYfullPred.Zs_non       = repmat(simple.z, 1, size(strat.Cs, 2));
    stationStruct.concProf.MYfullPred.Zs           = repmat(strat.Zs, 1, size(strat.Cs, 2));
    stationStruct.concProf.MYfullPred.Zs_nosed     = repmat(strat_nosed.Zs, 1, size(strat.Cs, 2));
    stationStruct.concProf.MYfullPred.ZsSum        = strat.ZsSum;
    
    stationStruct.concProf.MYfullPred.K_red = mean(stationStruct.velProf.MYfullPred.nu_t) / mean(stationStruct.velProf.MYfullPred.nu_t_nosed);

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
            plot(stationStruct.Velocity.data.velMean, stationStruct.Velocity.data.measZ./h, 'ko')
            xlabel('$u$ (m/s)')
            ylabel('$z/h$')
        subplot(1,4,2); hold on;
            cla
            plot(simple.c_non.*sum(cb_i), simple.eta, 'k--')
            plot(strat.c_hat.*strat.C_mean, strat.eta, 'r')
            plot(strat_nosed.c_hat.*strat_nosed.C_mean, strat_nosed.eta, 'b')
            plot(stationStruct.waterSamplesTable.concNW./2650, stationStruct.waterSamplesTable.sampleZ./h, 'ko')
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
