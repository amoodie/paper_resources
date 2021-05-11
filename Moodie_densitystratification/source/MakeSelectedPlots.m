function MakeSelectedPlots()
    
    source_path = genpath('.');
    utils_path = genpath('utils');
    data_path = genpath(fullfile('..', 'data'));
    addpath(source_path, utils_path, data_path);
    
    load('toSelectedPlots.mat')
    con = load_conset('quartz-water');

    if isempty(df(1).velProf)
        error('no velProf listing, are you sure you ran MakeStratModels since running MakeDataTable?')
    end

    global colorSet col1 col2 col3 col4 col5 mar1 mar2 mar3
    colswitch = 'colorful';
    [colorSet, col1, col2, col3, col4] = load_colorSet(colswitch);
    mar1 = 's';
    mar2 = 'o';
    mar3 = '^';
    markers = {mar1; mar2; mar3};  
    bgcolor = [1 1 1];
    col5 = [0.482 0.3059 0.5843];
    
    WP = readtable('wright_dissertation_combined.csv');
    WP.Rep = (sqrt(con.R .* con.g .* (WP.D50_mm_./1e3)) .* (WP.D50_mm_./1e3)) ./ con.nu;
    WP.Dstar = ((con.R * con.g) / (con.nu * con.nu))^(1/3) .* (WP.D50_mm_./1e3);
    WP.river(strcmp(WP.river,'atchafalaya')) = {'Atchafalaya'};
    WP.river(strcmp(WP.river,'midloupe')) = {'Mid Loup'};
    WP.river(strcmp(WP.river,'mississippi')) = {'Mississippi'};
    WP.river(strcmp(WP.river,'niobrara')) = {'Niobrara'};
    WP.river(strcmp(WP.river,'red')) = {'Red'};
    WP.river(strcmp(WP.river,'riogrande')) = {'Rio Grande'};
    
    waterSampleZs =  [0.05, 0.15, 0.25, 0.5, 0.9];
    sampleZs =  [0, 0.05, 0.15, 0.25, 0.5, 0.9];
    
    % mesh for ustar for compatibility
    load('stratificationFieldDataParameterSpace.mat', 'out', 'K_red', 'field_data')
    pt = array2table([field_data.Ent, field_data.S, field_data.D50_m, field_data.Hbf_m, ...
                  field_data.ustar, K_red(:), field_data.ws], ...
                 'VariableNames', {'cb', 'S', 'D', 'h', 'ustar', 'K_red', 'ws'});
    pt = pt(pt.D < 2000e-6, :);
    
    % make or load the viscosit profiles
    if false
        viscProfs.H = 5;
        viscProfs.ustar = 0.1;
        viscProfs.D_dim = 90e-6;
        viscProfs.cb_i = [0.001];
        viscProfs.S_0 = 6.4e-5;
        viscProfs.Zs = linspace(0.05*viscProfs.H, viscProfs.H, 51);
        viscProfs.etas = viscProfs.Zs ./ viscProfs.H;
        viscProfs.paraClear = 0.41 .* viscProfs.ustar .* viscProfs.Zs .* (1-viscProfs.etas);
        [caseStruct] = MellorYamadaConcVelModel(viscProfs.D_dim, viscProfs.cb_i, viscProfs.ustar, viscProfs.S_0, viscProfs.H);
        viscProfs.MYclear = caseStruct.nu_t_nosed;
        viscProfs.MYsed = caseStruct.nu_t;
        viscProfs.alpha = calculate_alpha_WP04(viscProfs.cb_i, viscProfs.S_0);
        viscProfs.paraAlpha = viscProfs.alpha .* 0.41 .* viscProfs.ustar .* viscProfs.Zs .* (1-viscProfs.etas);
        save('../data/viscositProfiles.mat', 'viscProfs')
    else
        load('viscositProfiles.mat', 'viscProfs')
    end
    
    % combt for regression
    combt = array2table([ [(1 ./ WP.vs50_u_); station.suspension.suspensionNumberBedDistWeighted], ...
                          [WP.Rep; station.suspension.bedDistWeightedRep], ...
                          [WP.So; station.slope], ...
                          [WP.C5t; station.cb5S.cb5], ...
                          [WP.Kred; station.alpha.byPredD50RouToGsClassModelRou], ...
                          [ones(size(WP.Kred)); zeros(size(station.slope))] ], 'VariableNames', ...
                          {'ustarws', 'Rep', 'S', 'cb', 'K_red', 'WP'});
    combt.logS = log10(combt.S);
    combt.logcb = log10(combt.cb);
    combt.WP_alpha_pred = compute_wright_parker(combt.cb, combt.S);
    
    
    % indicies vectors
    idx.WPOut = true(size(WP.river));
    idx.velProf = arrayfun(@(x) ~(size(fieldnames(x.velProf),1)<5), df);
    
    % predictor table fromm model data]
    pt.ustarws = pt.ustar ./ pt.ws;
    pt.Rep = (sqrt(con.g .* con.R .* (pt.D)) .* pt.D) ./ con.nu;
    pt.Dstar = ((con.R * con.g) / (con.nu * con.nu))^(1/3) .* pt.D;
    pt.logS = log10(pt.S);
    pt.logcb = log10(pt.cb);
    
    % printing options
    printOut = true;
    print_root = '../figures/';
    
    
    %% velocity calibration plot(s)    
    % main demonstration plot for the calibration, goes into paper
    ustarActList = [station.velCalib.ustarIntoCalib(station.velCalib.ustarYr==2016); station.velCalib.ustarIntoCalib(station.velCalib.ustarYr==2018)];
    ustarPredList = [station.velCalib.ustarOutofCalib(station.velCalib.ustarYr==2016); station.velCalib.ustarOutofCalib(station.velCalib.ustarYr==2018)];
    ustar_mdl_r2 = 1 - (sum((ustarActList - ustarPredList).^2)) / ((length(ustarActList)-1) * var(ustarActList));
    before_correlation = corrcoef((nearBed.FlowDepth), nearBed.concNW);
    after_correlation = corrcoef(nearBed.ustarCalib, nearBed.concNW);
    
    f = figure();
    subplot(1, 3, 1); hold on;
        [YRmarkers] = plotYRdata(nearBed.FlowDepth, nearBed.concNW, idx.sampleYr, idx.sampleOut, f);
        text(0.6, 0.1, ['r $=' num2str(before_correlation(2, 1), '%.2f'), '$'], 'Units', 'Normalized')
        xlabel('flow depth $H$ (m)')
        ylabel('near bed conc. $\bar{c}_b$ (g/L)')
        box on
        axis square
        legend(YRmarkers, {'2015', '2016', '2018'}, 'EdgeColor', [0 0 0], 'Interpreter', 'latex')
        text(0.05, 0.9, 'a', 'units', 'normalized', 'FontSize', 14)
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot(1, 3, 2); hold on;
        plot(station.velCalib.ustarIntoCalib(station.velCalib.ustarYr==2016), station.velCalib.ustarOutofCalib(station.velCalib.ustarYr==2016), ...
            'Marker', mar2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col2, 'MarkerSize', 6, 'LineStyle', 'none')
        plot(station.velCalib.ustarIntoCalib(station.velCalib.ustarYr==2018), station.velCalib.ustarOutofCalib(station.velCalib.ustarYr==2018), ...
            'Marker', mar3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col3, 'MarkerSize', 6, 'LineStyle', 'none')
        plot([0 0.1], [0, 0.1], 'k-')
        text(0.07, 0.08, '1:1 line', 'rotation', 45)
        text(0.6, 0.1, ['R$^2=' num2str(ustar_mdl_r2, '%.2f'), '$'], 'Units', 'Normalized')
        axis square
        ylabel('$u_*$ calibrated (m/s)')
        xlabel({'$u_*$ from friction relation (m/s)'})
        xlim([0 0.1])
        ylim([0 0.1])
        box on
        text(0.05, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot(1, 3, 3); hold on;
        plotYRdata(nearBed.ustarCalib, nearBed.concNW, idx.sampleYr, idx.sampleOut, f);
        text(0.6, 0.1, ['r $=' num2str(after_correlation(2, 1), '%.2f'), '$'], 'Units', 'Normalized')
        axis square
        ylabel('near bed conc. $\bar{c}_b$ (g/L)')
        xlabel('$u_*$ calibrated (m/s)')
        box on
        text(0.05, 0.9, 'c', 'units', 'normalized', 'FontSize', 14)
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    pause(0.5);
    set(gcf, 'Pos', [50 100 950 300], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    drawnow; pause(1);
    if printOut
        print('-dpng', '-r300', [print_root, 'shear_velocity_calibration.png']);
        set(gcf, 'InvertHardcopy', 'on')
        print('-depsc', '-painters', '-r300', [print_root, 'shear_velocity_calibration.eps']);
        pause(1);
    end
    
    % a supplementary plot to demonstrate how the new u* differs from the DSP
    f = figure();
        hold on;
        [mkData] = plotYRdata(station.depth, station.ustarCalib, idx.stationYr, idx.stationOut, f);
        [stdDSP] = plot(0:0.2:8, (sqrt(con.g .* (0:0.2:8) .* df(1).Velocity.slope)), 'k--');
        xlabel('depth (m)')
        ylabel('$u_*$ calibrated (m/s)')
        box on
        legend([stdDSP], {'depth-slope product'}, 'Location', 'SouthEast', 'EdgeColor', [0 0 0], 'Interpreter', 'latex')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 300 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        if printOut
            print('-dpng', '-r300', [print_root, 'depthVsUstarCalibration.png']);
            set(gcf, 'InvertHardcopy', 'on')
            print('-depsc', '-painters', '-r300', [print_root, 'depthVsUstarCalibration.eps']);
            pause(1);
        end
      
        
    %% viscosity profiles demo
    f = figure(); hold on;
        %[clear] = plot(viscProfs.paraClear, viscProfs.etas, '--', 'Color', col4, 'MarkerSize', 3, 'LineWidth', 2);
        %[wpalpha] = plot(viscProfs.paraAlpha, viscProfs.etas, '-', 'Color', col4, 'MarkerSize', 3, 'LineWidth', 2);
        [MYsum_mdl] = plot(viscProfs.MYsed, viscProfs.etas, '-', 'Color', col5, 'LineWidth', 1.5); % mellor yamada prediction with sed
        [MYsum_mdl_clear] = plot(viscProfs.MYclear, viscProfs.etas, '--', 'Color', col5, 'LineWidth', 1.5); % mellor yamada prediction clear
        xlabel('eddy viscosity')
        ylabel('$z/H$')
        box on
        legend([MYsum_mdl_clear, MYsum_mdl], {'$MY$ no sediment', '$MY$'},...
           'Location', 'NorthEast', 'EdgeColor', [0 0 0], 'interpreter', 'latex')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 400 400], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        if printOut
            print('-dpng', '-r300', [print_root, 'viscosityProfiles.png']);
            set(gcf, 'InvertHardcopy', 'on')
            print('-depsc', '-painters', '-r300', [print_root, 'viscosityProfiles.eps']);
            pause(1);
        end
        
    f = figure; hold on;
        plot([0 0.15], [0, 0.15], 'k-')    
        [YRmarkers] = plotYRdata(station.velCalib.ustarFric, station.velCalib.ustarFric_MaSolved, ...
            idx.sampleYr(idx.velProf), idx.sampleOut(idx.velProf), f);
        xlabel('$u_*$ from friction relation (Moodie solved) (m/s)')
        ylabel('$u_*$ from friction relation (Ma solved) (m/s)')
        text(0.1, 0.12, '1:1 line', 'rotation', 45)
        box on
        axis square
        legend(YRmarkers, {'2016', '2018'}, 'EdgeColor', [0 0 0], 'Interpreter', 'latex', 'location', 'southeast')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')



    %% theta and alpha plots
    nClasses = size(df(1).nearBedData.gsDistNearBedNWnorm, 1);   
    f = figure();
    subplot(2, 1, 1); hold on; % theta plot for near bed distribution ws
        [YRmarkers] = plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted, station.theta.meanNormToTotalModel, idx.stationYr, idx.stationOut, f);
        plot([1e-1 1e4], [0 0], 'k:', 'LineWidth', 1.5)
        text(0.05, 0.9, 'a', 'units', 'normalized', 'FontSize', 14)
        xlabel('dimensionless shear velocity $(u_* / w_{s})$')
        ylabel('$\hat{\theta} / \bar{c}_b$')
        xlim([0 25])
        ylim([-0.5 0.5])
        box on
        legend(YRmarkers, {'2015', '2016', '2018'}, 'Location', 'SouthEast', 'EdgeColor', [0 0 0], 'Interpreter', 'latex')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot(2, 1, 2); hold on;
        [markerset2] = plotYRdata(station.suspension.suspensionNumberGsClass, station.theta.meanNormGsClassToGsClassModel, ...
            varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, 'ColorByCategorical', station.gsClassVectorCats);
        plot([1e-1 1e4], [0 0], 'k:', 'LineWidth', 1.5)
        text(0.05, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
        set(gca, 'xscale', 'log')
        xlabel('dimensionless shear velocity ($u_* / w_{s,i}$)')
        ylabel('$\hat{\theta}_i / \bar{c}_{b,i}$')
        colb = colorbar;
        cleanup_colorbar(colb, round(station.gsClassVector(1:nClasses),0), ...
            'title', 'grain-size class ($\mu$m)')
        colormap((parula(nClasses)))
        xlim([5e-1 5e2])
        ylim([-1 1])
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        pause(1);
        if printOut
            set(gcf, 'Pos', [50 100 600 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
            print('-dpng', '-r300', [print_root, 'normalized_mean_signed_deviation.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'normalized_mean_signed_deviation.eps']);
            pause(1);
        end   
                
    f = figure(); hold on; % alpha with near bed suspension with Rep
        plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted .* station.suspension.nearBedNWnormDistWeightedRep .^ 0.6, station.alpha.byPredD50RouToGsClassModelRou, ...
            idx.stationYr, idx.stationOut, f);
        plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
        ylim([0 1.5])
        xlabel('$(u_* / w_{s,D_{w},nb})$ Rep$^{0.6}$')
        ylabel('$\alpha_f$')
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        if printOut
            print('-dpng', '-r300', [print_root, 'suspensionRep_YRonly.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'suspensionRep_YRonly.eps']);
            pause(1);
        end
        
    f = figure(); % alpha in WP space from comparing rouse profiles
    subplot(2, 1, 1); hold on;
        [cld] = fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
        [YRmarkers] = plotYRdata(station.cb5S.byGsClassSummedModel, station.alpha.byPredD50RouToTotalModelRou, idx.stationYr, idx.stationOut, f);
        plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
        plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
        text(0.05, 0.9, 'a', 'units', 'normalized', 'FontSize', 14)
        xlim([0 250])
        ylim([0 1.5])
        xlabel('$\bar{c}_{b} / S_0$')
        ylabel('$\alpha_{f}$')
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        legend([cld; YRmarkers], {'WP04 data extent', '2015', '2016', '2018'}, 'Location', 'NorthEast', 'EdgeColor', [0 0 0], 'Interpreter', 'latex')
    subplot(2, 1, 2); hold on;
        fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
        plotYRdata(station.cb5S.byGsClassModelToGsClassFrac, station.alpha.byPredGsClassRouToGsClassModelRou, ...
            varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, 'ColorByCategorical', station.gsClassVectorCats);
        plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
        colb = colorbar;
        cleanup_colorbar(colb, round(station.gsClassVector(1:nClasses),0), ...
            'title', 'grain-size class ($\mu$m)')
        colormap((parula(nClasses)))
        plot([0 300], [1 1], 'k:', 'LineWidth', 1.5)
        text(0.05, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
        set(gca, 'yscale', 'log')
        ylim([1e-1, 5e1])
        xlim([0 250])
        xlabel('$(\bar{c}_{b,i}/F_i) / S_0$')
        ylabel('$\alpha_{f,i}$')
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        pause(0.5)
        set(gcf, 'Pos', [50 100 600 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        if printOut
            print('-dpng', '-r300', [print_root, 'alpha_data_plots.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'alpha_data_plots.eps']);
            pause(1);
        end
        
    f = figure();
    subplot(2, 1, 1); hold on; % comparison of alha from field data to alpha from MY predictions
        plotYRdata(station.alpha.byPredD50RouToGsClassModelRou, station.alpha.byWP04func, idx.stationYr, idx.stationOut, f);
        plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
        plot([0 1], [0 1], 'k-', 'LineWidth', 1)
        text(0.05, 0.9, 'a', 'units', 'normalized', 'FontSize', 14)
        ylim([0 1])
        xlim([0 1])
        axis square
        xlabel('$\alpha_f$')
        ylabel('$\alpha \equiv \alpha_{WP04}$')
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot(2, 1, 2); hold on; % comparison of alha from field data to alpha from MY predictions
        plotYRdata(station.alpha.byPredD50RouToGsClassModelRou, station.alpha.byKredMY, idx.stationYr, idx.stationOut, f);
        plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
        plot([0 1], [0 1], 'k-', 'LineWidth', 1)
        text(0.05, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
        ylim([0 1])
        xlim([0 1])
        axis square
        xlabel('$\alpha_f$')
        ylabel('$\alpha \equiv K_{red,MY}$')
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 600 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        if printOut
            print('-dpng', '-r300', [print_root, 'alpha_data_correlations.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'alpha_data_correlations.eps']);
            pause(1);
        end
        
    f = figure(); hold on;
    subplot(1, 2, 1); hold on;
        YRmarkers = plotYRdata(station.alpha.byPredD50RouToGsClassModelRou, station.Risr.byNearBedBulk, idx.stationYr, idx.stationOut, f);
        WPmarkers = plotWPdata(WP.Kred, WP.Rib, categorical(WP.river), idx.WPOut, f);
        text(0.9, 0.9, 'a', 'units', 'normalized', 'FontSize', 14)
        xlabel('$\alpha_f$ or $\overline{K}_{red}$')
        ylabel('sand-river Richardson num. $\textrm{Ri}_{\textrm{sr}}$')
        plot([1 1], [1e-2, 1e2], 'k:', 'LineWidth', 1.5)
        ylim([1e-2, 3e1])
        xlim([0, 1.25])
        box on
        axis square
        lgd = legend([YRmarkers; WPmarkers], vertcat({'2015'}, {'2016'}, {'2018'}, unique(WP.river)), 'location', 'southwest', 'Interpreter', 'latex');
        set(gca, 'yscale', 'log')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot(1, 2, 2); hold on;
        plot([1 1], [1 100], 'k:', 'LineWidth', 1.5)
        plotYRdata(station.alpha.byPredGsClassRouToGsClassModelRou, station.beta.fieldr0Long, ...
            varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, ...
            'ColorByCategorical', station.gsClassVectorCats)
        % below would plot the nearbed concentration on the xaxis
        %plotYRdata(station.cb5S.byGsClassModelToGsClassFrac.*6.4e-5, station.beta.fieldr0Long, ...
        %    varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, ...
        %    'ColorByCategorical', station.gsClassVectorCats)
        text(0.9, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
        colb = colorbar;
        cleanup_colorbar(colb, round(station.gsClassVector(1:nClasses),0), ...
            'title', 'grain-size class ($\mu$m)')
        colormap((parula(nClasses)))
        set(gca, 'xscale', 'log', 'yscale', 'log')
        xlabel('$\alpha_{f,i}$')
        ylabel('recovery coeff. $r_{0,i}$')
        ylim([1 100])
        xlim([5e-2 5e1])
        xticks([1e-1, 1e0])
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 800 350], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        if printOut
            print('-dpng', '-r300', [print_root, 'alpha_vs_sandrichardson.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'alpha_vs_sandrichardson.eps']);
            pause(1);
        end

    
    %% sediment diffusivity alphas and ratio plot
    m_plot = 0.1; % margins
    s_plot = 0.125; % spacing
    nx = 2;
    ny = 2;
    w_plot = (1-2*m_plot-(nx-1)*s_plot) / nx;
    h_plot = (1-2*m_plot-(ny-1)*s_plot) / ny;
    x_pos = repmat(m_plot+(0:(nx-1))*w_plot+(0:(nx-1))*s_plot, ny, 1)';
    y_pos = repmat(1-(h_plot+m_plot+(0:(ny-1))*h_plot+(0:(ny-1))*s_plot), nx, 1);
    f = figure(); hold on;
    sp1 = subplot('position', [x_pos(1) y_pos(1) w_plot h_plot]); hold on;
        plot([0 2], [0 2], 'k-', 'LineWidth', 1);
        plot([1 1], [0 2], 'k:', 'LineWidth', 1.5);
        plot([0 2], [1 1], 'k:', 'LineWidth', 1.5);
        [YrMarks] = plotYRdata(station.alpha.byKredVelocityModel, station.alpha.byPredD50RouToGsClassModelRou(idx.velProf), ...
            idx.stationYr(idx.velProf), idx.stationOut(idx.velProf), f);
        text(0.05, 0.9, 'a', 'units', 'normalized', 'FontSize', 14)
        xlabel('$\alpha_{V}$')
        ylabel('$\alpha_f$')
        xlim([0 1.5]); ylim([0 1.5]);
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    sp2 = subplot('position', [x_pos(3) y_pos(3) w_plot h_plot]); hold on;
        plot([1e-1 200], [1 1], 'k:', 'LineWidth', 1.5);
        plotYRdata(station.beta.byGsClassModelToGsClassFrac, station.beta.alphaGsClassVelProfOnly .* (1./station.beta.alphaByVelocityVelProfOnly), ...
            repelem(idx.stationYr(idx.velProf), nClasses), repelem(idx.stationOut(idx.velProf), nClasses), f, ...
            'ColorByCategorical', station.gsClassVectorCats(1:(sum(idx.velProf)*nClasses)));
        text(0.875, 0.875, 'b', 'units', 'normalized', 'FontSize', 14);
        colb = colorbar;
        cleanup_colorbar(colb, round(station.gsClassVector(1:nClasses),0), ...
            'title', 'grain-size class ($\mu$m)')
        colormap((parula(nClasses)))
        xlabel('$\bar{c}_{b,i}/F_i$')
        ylabel('$\beta_i$')
        xlim([1e0 100])
        set(gca, 'xscale', 'log', 'yscale', 'log')
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot('position', [x_pos(4) y_pos(4) w_plot 2*h_plot+s_plot]); hold on
        plot([1 1], [0 1], '--', 'Color', 'k', 'LineWidth', 1);
        plot([0 0], [0 1], ':', 'Color', 'k', 'LineWidth', 1.5);
        [YrMarks] = plotYRdata(sdif.tab.tab.pe, sdif.tab.tab.zh+(0.06*(rand(size(sdif.tab.tab.zh))-0.5)), ...
                sdif.tab.tab.cyear, ones(size(sdif.tab.tab.cyear)), f, 'ColorByCategorical', ...
                categorical(cellstr(num2str(round(sdif.tab.tab.gs)))), 'ColorEdgeOnly', true);
        text(0.9, 0.925, 'c', 'units', 'normalized', 'FontSize', 14)
        colb = colorbar;
        cleanup_colorbar(colb, round(station.gsClassVector(1:nClasses),0), ...
            'title', 'grain-size class ($\mu$m)');
        colormap((parula(nClasses)))
        xlabel('$(\bar{c}_{f,i} - c_{MY,i})/\bar{c}_{f,i}$')
        ylabel('$z/H$')
        xlim([-2 2])
        set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.50','0.9'})
        %axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 650 450], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        pause(0.5); set(sp1, 'Position', [sp2.Position(1), y_pos(1), sp2.Position(3), h_plot])
        if printOut
            print('-dpng', '-r300', [print_root, 'sediment_diff.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'sediment_diff.eps']);
            pause(1);
        end
    
    % supplementary plot for beta by grain size
    f = figure(); hold on 
        plot([1e-1 200], [1 1], 'k:', 'LineWidth', 1.5)
        plotYRdata(station.beta.byGsClassModelToGsClassFrac, station.beta.alphaGsClassVelProfOnly .* (1./station.beta.alphaByVelocityVelProfOnly), ...
            repelem(idx.stationYr(idx.velProf), nClasses), repelem(idx.stationOut(idx.velProf), nClasses), f, ...
            'ColorByCategorical', station.gsClassVectorCats(1:(sum(idx.velProf)*nClasses)))
        colb = colorbar;
        cleanup_colorbar(colb, round(station.gsClassVector(1:nClasses),0), ...
            'title', 'grain-size class ($\mu$m)')
        colormap((parula(nClasses)))
        xlabel('$\bar{c}_{b,i}/F_i$')
        ylabel('$\beta_i$')
        xlim([1e0 100])
        set(gca, 'xscale', 'log', 'yscale', 'log')
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 400 350], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        if printOut
            print('-dpng', '-r300', [print_root, 'beta_grainsize.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'beta_grainsize.eps']);
            pause(1);
        end

        % do the cumulate as a supplementary figure (to have both isolated)
        f = figure(); hold on;
            plot([0 50], [1 1], 'k:', 'LineWidth', 1.5)
            plotYRdata(station.cb5S.cb5(idx.velProf).*2650, station.alpha.byPredD50RouToGsClassModelRou(idx.velProf).*(1./station.alpha.byKredVelocityModel), ...
                idx.stationYr(idx.velProf), idx.stationOut(idx.velProf), f)
            text(0.875, 0.875, 'b', 'units', 'normalized', 'FontSize', 14)
            xlabel('$\bar{c}_b$ (g/L)')
            ylabel('$\beta$')
            xlim([0 50]); ylim([0 1.5]);
            axis square
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 400 350], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        if printOut
            print('-dpng', '-r300', [print_root, 'beta_cumulated.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'beta_cumulated.eps']);
            pause(1);
        end
    
    %% alpha with gamma correction plots -- WPalpha - fit data
    cmap = parula(nClasses); %flipud(jet(6));
    m_plot = 0.1; % margins
    s_plot = 0.125; % spacing
    nx = 2;
    ny = 2;
    w_plot = (1-2*m_plot-(nx-1)*s_plot) / nx;
    h_plot = (1-2*m_plot-(ny-1)*s_plot) / ny;
    x_pos = repmat(m_plot+(0:(nx-1))*w_plot+(0:(nx-1))*s_plot, ny, 1)';
    y_pos = repmat(1-(h_plot+m_plot+(0:(ny-1))*h_plot+(0:(ny-1))*s_plot), nx, 1);
    figure();
    subplot(2, 2, [1 3]); hold on
        plot([0 0], [0 max(reg.Zs)], '-', 'Color', [0.6 0.6 0.6])
        sdd_h = plot(reg.sdd(1:49,:).*2650, reg.Zs(1:49), 'LineWidth', 2);
        set(sdd_h, {'color'}, num2cell(cmap, 2));
        text(0.1, 0.95, 'a', 'units', 'normalized', 'FontSize', 14)
        ylim([0 max(reg.Zs)])
        xlim([-0.5 1])
        ylabel('$z/H$ (-)')
        xlabel('$\xi_{WP04,i}$ (g/L)')
        cleanup_boxplot('color', jet(0))
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot('position', [x_pos(2) y_pos(2) w_plot, h_plot]); hold on
        plot([0 7], [0 0], '-', 'Color', [0.6 0.6 0.6])
        text(0.9, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
        boxplot(squeeze(reg.msd)'.*2650)
        cleanup_boxplot('color', cmap)
        set(gca, 'xticklabels', arrayfun(@(x) num2str(round(x,0)), gsDistClass2num(df(1).nearBedData.gsDistNearBedNWnorm), 'unif', 0))
        xlabel('grain-size ($\mu$m)')
        ylabel('$\int_b^H \xi_{WP04,i}/H$ (g/L)')
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot('position', [x_pos(4) y_pos(4) w_plot, h_plot]); hold on
        plot([0 7], [0 0], '-', 'Color', [0.6 0.6 0.6])
        text(0.9, 0.9, 'c', 'units', 'normalized', 'FontSize', 14)
        boxplot(squeeze(reg.nmsd)')
        cleanup_boxplot('color', cmap)
        set(gca, 'xticklabels', arrayfun(@(x) num2str(round(x,0)), gsDistClass2num(df(1).nearBedData.gsDistNearBedNWnorm), 'unif', 0))
        xlabel('grain-size ($\mu$m)')
        ylabel('$[\int_b^H \xi_{WP04,i}/H]/\bar{c}_b$ (-)')
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 600 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    if printOut
        print('-dpng', '-r300', [print_root, 'WPalpha_grainsize_msd.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'WPalpha_grainsize_msd.eps']);
        pause(1);
    end
    
    %% alpha regression plots
    % alpha regression for slopes
    Rep_exp = -0.6;
    
    x0 = [1, 52, 10, 31];
    min_func = (  @(x) minimize_regress(x, combt.K_red, combt.cb .* combt.Rep .^ Rep_exp, combt.S)  );
  
    xp = fminsearch(min_func, x0);
    w = xp(1);
    x = xp(2);
    y = xp(3);
    z = xp(4);
       
    rndw = round(w, 2, 'significant'); 
    rndx = round(x, 2, 'significant'); 
    rndy = round(y, 2, 'significant'); 
    rndz = round(z, 2, 'significant');
 
    mfunc = @(S) rndw .* log10(S) + rndx;
    nfunc = @(S) rndy .* log10(S) + rndz;
    Yfunc = @(m, n, x) 1 - ( m .* x .^ n);
    
    Ss = logspace(-5, -3, 5);
    Xs = logspace(-8, 0, 100)';
    Ys = Yfunc(mfunc(Ss), nfunc(Ss), Xs);
    
    XData1 = (pt.cb) .* pt.Rep.^Rep_exp;
    XData2 = (station.cb5S.cb5) .* station.suspension.bedDistWeightedRep.^Rep_exp;
    XData3 = (WP.C5t) .* WP.Rep.^Rep_exp;
    
    % total lists
    listPred = [Yfunc(mfunc(pt.S), nfunc(pt.S), XData1); Yfunc(mfunc(station.slope), nfunc(station.slope), XData2);  Yfunc(mfunc(WP.So), nfunc(WP.So), XData3)];
    listAct = [pt.K_red; station.alpha.byPredD50RouToGsClassModelRou; WP.Kred];
    slope_mdl_r2 = 1 - (nansum((listAct - listPred).^2)) / ((length(listAct)-1) * var(listAct, 'omitnan'));
    
    f = figure(); 
    sp1 = subplot(1,2,1); hold on;
        plot(Xs, Ys, 'k-', 'MarkerSize', 4)
        scatter(pt.cb.*pt.Rep.^Rep_exp, pt.K_red, 55, pt.S, 'filled', 'MarkerFaceAlpha', 0.5);
        set(gca, 'XScale', 'log')
        plotYRdata(station.cb5S.cb5.*station.suspension.bedDistWeightedRep.^Rep_exp, station.alpha.byPredD50RouToGsClassModelRou, idx.stationYr, idx.stationOut, f, ...
            'ColorByVariable', station.cb5S.S);
        [WP_marks] = plotWPdata(WP.C5t.*WP.Rep.^Rep_exp, WP.Kred, categorical(WP.river), idx.WPOut, f, ...
            'ColorByVariable', WP.So);
        text(Xs(57), 0.5, sprintsci(Ss(1),2,'Interpreter','latex'), 'rotation', -63)
        text(Xs(68), 0.54, sprintsci(Ss(3),2,'Interpreter','latex'), 'rotation', -56)
        text(Xs(65), 0.81, sprintsci(Ss(5),2,'Interpreter','latex'), 'rotation', -28)
        text(0.05, 0.5, 'a', 'units', 'normalized', 'FontSize', 14)
        ylim([0, 1])
        xlim([1e-6 1e-2])
        xlabel(['$\bar{c}_b \textrm{Re}_{\textrm{p,50}}^{', num2str(Rep_exp), '}$'])
        ylabel('$\alpha_{S_0}$')
        cb = colorbar;
        colormap(plasma())
        cb.Label.String = '${S_0}$'; cb.Label.Interpreter = 'latex';
        caxis([1e-5 1e-3])
        axis square
        set(gca, 'colorscale', 'log')
        box on
        legend(WP_marks, unique(WP.river), 'location', 'southwest', 'Interpreter', 'latex')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    sp2 = subplot(1, 2, 2); hold on;
        scatter(pt.K_red, Yfunc(mfunc(pt.S), nfunc(pt.S), XData1), 55, pt.S, 'filled', 'MarkerFaceAlpha', 0.5);
        plotYRdata(station.alpha.byPredD50RouToGsClassModelRou, Yfunc(mfunc(station.slope), nfunc(station.slope), XData2), ...
            idx.stationYr, idx.stationOut, f, 'ColorByVariable', station.cb5S.S);
        [WP_marks] = plotWPdata(WP.Kred, Yfunc(mfunc(WP.So), nfunc(WP.So), XData3), categorical(WP.river), idx.WPOut, f, ...
            'ColorByVariable', WP.So);
        text(0.7, 0.1, ['R$^2=' num2str(round(slope_mdl_r2, 2)), '$'])
        text(0.05, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
        colormap(plasma())
        caxis([1e-5 1e-3])
        set(gca, 'colorscale', 'log')
        plot([0 1], [0 1], 'k-')
        xlim([0 1]);
        ylim([0 1]);
        xlabel('$\alpha_f$ or $\overline{K}_{red}$')
        ylabel('$\alpha_{S_0}$ predicted')
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    set(gcf, 'Pos', [50 100 1000 375], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor);pause(0.5);
    sp2.Position = [sp1.Position(1) + 0.5, sp1.Position(2:4)];
    if printOut
        print('-dpng', '-r300', [print_root, 'alpha_regression3.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'alpha_regression3.eps']);
        pause(1);
    end
    
    
    % alpha regression for grain size
    regTabTot = array2table([(station.suspension.suspensionNumberGsClass), (station.cb5S.byGsClassModelToGsClassFrac), ...
        (station.alpha.byPredGsClassRouToGsClassModelRou), log(station.suspension.suspensionNumberGsClass), ...
        log(station.cb5S.byGsClassModelToGsClassFrac), log(station.alpha.byPredGsClassRouToGsClassModelRou), ...
        log((station.suspension.suspensionNumberGsClass) .* station.cb5S.byGsClassModelToGsClassFrac), ...
        station.gsClassVector, varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses)], ...
        'VariableNames', {'ustarws', 'cb5S', 'alpha', 'logustarws', 'logcb5S', 'logalpha', 'logustarwscb5S', 'D', 'idxYr', 'idxOut'});
    regTab = regTabTot;
    regTab(or(isinf(regTab.logalpha), isinf(regTab.logcb5S)), :) = [];
    regTab(~regTab.idxOut, :) = [];
    regTabTot.idxOut(or(isinf(regTab.logalpha),isinf(regTab.logcb5S))) = 0;
    regmdl = fitnlm(regTab, 'logalpha ~ a * logustarwscb5S + b', [-0.3  0.4]);
    arnd = round(exp(regmdl.Coefficients.Estimate(2)),2, 'significant');
    brnd = round(regmdl.Coefficients.Estimate(1), 2, 'significant');
    reg_mdl_r2 = regmdl.Rsquared.Ordinary;
    
    f = figure(); hold on;
        plot([1e-4, 100], [1 1], 'k:', 'LineWidth', 1.5)
        plot([1e-4, 100], [1e-4, 100], 'k-')
        plotYRdata(regTabTot.alpha, (regTabTot.ustarws .* regTabTot.cb5S).^brnd .* arnd, ...
            regTabTot.idxYr, regTabTot.idxOut, f, 'ColorByCategorical', station.gsClassVectorCats);
        set(gca, 'xscale', 'log', 'yscale', 'log')
        xlim([1e-2, 100])
        ylim([1e-2, 100])
        ylabel(['$\alpha_{\textrm{YR}} \equiv m \left[\frac{u_*}{w_s}\frac{\bar{c}_b}{S_0}\right]^n$'])
        xlabel('$\alpha_{f,i}$')
        text(0.7, 0.1, ['R$^2=' num2str(round(reg_mdl_r2, 2)), '$'], 'Units', 'Normalized')
        axis square
        colb = colorbar;
        colormap((parula(nClasses)))
        cleanup_colorbar(colb, round(station.gsClassVector(1:nClasses),0), ...
            'title', 'grain-size class ($\mu$m)')
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 500 300], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    if printOut
        print('-dpng', '-r300', [print_root, 'alpha_grainsize_regression.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'alpha_grainsize_regression.eps']);
        pause(1);
    end
   
    %% entrainment plot
    % aux calculations for text on entrainment
    eAlpha1 = concatenate_cell_vectors(varRepLong(arrayfun(@(x) repmat(x, 1, nClasses-1), station.alpha.byPredD50RouToTotalModelRou, 'Unif', 0), nearBed.dims));
    w = entr.dist.EsiPredLongCourseOnly ./ entr.dist.EsiMeasLongCourseOnly;
    high_idx = w > 2;
    low_idx = w < 0.5;
    in_idx = ~or(high_idx, low_idx);
    cats = cell(size(high_idx));
    cats(high_idx) = {'h'}; cats(low_idx) = {'l'}; cats(in_idx) = {'i'};
    entr_grainsize = concatenate_cell_vectors(varRepLong(arrayfun(@(x) gsDistClass2num(x.bedData.gsDistBed), df, 'unif', 0), nearBed.dims));
    
    f = figure();
    subplot(1,2,1)
    hold on; set(gca, 'XScale', 'log', 'YScale', 'log')
        plotYRdata(entr.dist.lXiLong, entr.dist.EsiMeasLong, repelem(idx.sampleYr, 6, 1), repelem(idx.sampleOut, 6, 1), f, ...
            'ColorByCategorical', categorical(entr_grainsize));
        lineWP04 = loglog(entr.mods.WP04Zu, entr.mods.WP04Es, 'Color', [0 0 0], 'LineWidth', 2);
        text(0.05, 0.9, 'a', 'units', 'normalized', 'FontSize', 14)
        ylim([1e-5 1e0])
        xlim([1e0 5e1])
        yticks([1e-4 1e-3 1e-2 1e-1 1e0])
        colb = colorbar;
        colormap(gca, (parula(nClasses)));
        cleanup_colorbar(colb, round(station.gsClassVector(1:nClasses),0), ...
            'title', 'grain-size class ($\mu$m)');
        xlabel('$\lambda X_i$')
        ylabel('$E_{s,i} \equiv \bar{c}_{b,i} / F_i$')
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    sp2 = subplot(1, 2, 2);
    hold on; set(gca, 'XScale', 'log', 'YScale', 'log')
        plot([1e-5, 100], [1e-5, 100], 'k-')
        plot([1e-5, 100], [1e-6, 10], 'k--')
        plot([1e-6, 100], [1e-5, 1000], 'k--')
        plotYRdata(entr.dist.EsiMeasLong, entr.dist.EsiPredLong, repelem(idx.sampleYr, 6, 1), repelem(idx.sampleOut, 6, 1), f, ...
            'ColorByVariable', concatenate_cell_vectors(varRepLong(arrayfun(@(x) repmat(x, 1, nClasses), ...
            station.alpha.byPredD50RouToTotalModelRou, 'Unif', 0), nearBed.dims)), 'IgnoreOutliers', true);
        text(0.05, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
        ylim([1e-6 1e0])
        xlim([1e-5 1e1])
        xticks([1e-4 1e-2 1e0])
        colb = colorbar;
        caxis([0 1]);
        colormap(gca, pink);
        colb.Label.Interpreter = 'latex';
        colb.Label.String = '$\alpha_{f}$';
        xlabel('measured $E_{s,i}$')
        ylabel('predicted $E_{s,i}$')
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    set(gcf, 'Pos', [50 100 800 300], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor); pause(0.5)
    sp2pos = get(sp2, 'position');
    subax1 = axes('Position', [sp2pos(1)+sp2pos(3)-0.11, sp2pos(2)+0.0985, 0.1, 0.3]); hold on;
        plot(subax1, [0 4], [1 1], 'k:')    
        boxplot(eAlpha1, categorical(cats))
        ylim([0, 1.5])
        cleanup_boxplot()
        ylabel(subax1, '$\alpha_f$')
        xlabel(subax1, 'meas./pred.')
        set(subax1, 'YTickLabels', {'0', '', '1', ''})
        set(subax1, 'XTickLabels', {'$<0.5$', '', '$>2$'})
        set(subax1, 'XAxisLocation','top')
        set(subax1, 'FontSize', 8)
    drawnow; pause(0.1);
    if printOut
        print('-dpng', '-r300', [print_root, 'entrainment.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'entrainment.eps']);
        pause(1);
    end
    
    
    f = figure(); hold on;
    plotYRdata(entr.dist.EsMeas, entr.dist.EsPred, idx.sampleYr, idx.sampleOut, f, ...
            'ColorByVariable', varRepLong(station.alpha.byPredD50RouToTotalModelRou, nearBed.dims), ...
            'IgnoreOutliers', true);
    
    
    
    %% sediment transport data
    f = figure(); hold on;
        plot([1 1], [0 3], 'k:', 'LineWidth', 1.5)
        plot([0 3], [1 1], 'k:', 'LineWidth', 1.5)
        plotYRdata(station.alpha.byPredD50RouToTotalModelRou(idx.velProf), station.transport.fieldNoStrat_ratio, ...
            idx.stationYr(idx.velProf), idx.stationOut(idx.velProf), f);
        xlim([0 1.5])
        ylim([0 2])
        xlabel('$\alpha_f$')
        ylabel('$q_s / q_{s,1.0}$')
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [50 100 620 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    drawnow; pause(0.1);
    if printOut
        print('-dpng', '-r300', [print_root, 'sedtransport.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'sedtransport.eps']);
        pause(1);
    end
        
    f = figure();
    subplot(2, 2, 1); hold on;
        plot([1 1], [0 3], 'k:', 'LineWidth', 1.5)
        plot([0 3], [1 1], 'k:', 'LineWidth', 1.5)
        [YRmarkers] = plotYRdata(station.alpha.byPredD50RouToTotalModelRou(idx.velProf), station.transport.fieldNoStrat_ratio_upper, ...
            idx.stationYr(idx.velProf), idx.stationOut(idx.velProf), f);
        lgd = legend(YRmarkers, vertcat({'2016'}, {'2018'}), 'location', 'southwest', 'Interpreter', 'latex'); 
        xlim([0 1.5])
        ylim([0 2])
        text(0.9, 0.9, 'a', 'units', 'normalized', 'FontSize', 14)
        xlabel('$\alpha_f$')
        ylabel({'$q_s / q_{s,1.0}$ in upper 20\%', 'of water column'})
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot(2, 2, 3); hold on;
        plot([1 1], [0 3], 'k:', 'LineWidth', 1.5)
        plot([0 3], [1 1], 'k:', 'LineWidth', 1.5)
        plotYRdata(station.alpha.byPredD50RouToTotalModelRou(idx.velProf), station.transport.concRatio(idx.velProf), ...
            idx.stationYr(idx.velProf), idx.stationOut(idx.velProf), f);
        xlim([0 1.5])
        ylim([0 2])
        text(0.9, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
        xlabel('$\alpha_f$')
        ylabel({'$\bar{c} / \bar{c}_{1.0}$ in upper 20\%', 'of water column'})
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot(2, 2, 4); hold on;
        plot([1 1], [0 2], 'k:', 'LineWidth', 1.5)
        plot([1e-3 1e2], [1 1], 'k:', 'LineWidth', 1.5)
        plotYRdata(station.alpha.byPredGsClassRouToGsClassModelRou(varRepLong(idx.velProf, nClasses)), ...
            station.transport.concRatioGsClass(varRepLong(idx.velProf, nClasses)), ...
            repelem(idx.stationYr(idx.velProf), nClasses), repelem(idx.stationOut(idx.velProf), nClasses), f, ...
            'ColorByCategorical', station.gsClassVectorCats(1:(sum(idx.velProf)*nClasses)));
        xlabel('$\alpha_f$')
        ylabel({'$\bar{c}_i / \bar{c}_{i,1.0}$ in upper 20\%', 'of water column'})
        set(gca, 'xscale', 'log')
        text(0.9, 0.9, 'c', 'units', 'normalized', 'FontSize', 14)
        colb = colorbar;
        cleanup_colorbar(colb, round(station.gsClassVector(1:nClasses),0), ...
            'title', 'grain-size class ($\mu$m)');
        colormap((parula(nClasses)))
        ylim([0, 2])
        xlim([1e-1, 1e1])
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    set(gcf, 'Pos', [50 100 890 660], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    drawnow; pause(0.1);
    if printOut
        print('-dpng', '-r300', [print_root, 'sedtransport_upper20_supp.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'sedtransport_upper20_supp.eps']);
        pause(1);
    end
    
    %% WrightParker original plots
    WP.PRED = compute_wright_parker(WP.C5t, WP.So);
    WP_original_r2 = 1 - (nansum((WP.Kred - WP.PRED).^2)) / ((length(WP.Kred)-1) * var(WP.Kred, 'omitnan'));
    combt_r2 = 1 - (nansum((combt.K_red - combt.WP_alpha_pred).^2)) / ((length(combt.K_red)-1) * var(combt.K_red, 'omitnan'));
    
    
    
    fivepctexceed = array2table([[3.5e-5, 5e-5, 7.8e-5, 8e-4, 1.1e-3, 1.9e-3]', ...
                                 [2.2, 2.7, 3.15, 0.8, 0.62, 0.2]', ...
                                 [0.265, 0.24, 0.22, 0.4, 0.42, 0.32]', ...
                                 [0.64, 0.63, 0.7, 0.8, 0.84, 0.94]'], ...
                                 'VariableNames', {'slope', 'cms0', 'uws', 'kred'});
    m_plot = 0.075; % margins
    s_plot = 0.1; % spacing
    nx = 2;
    ny = 2;
    w_plot = (1-2*m_plot-(nx-1)*s_plot) / nx;
    h_plot = (1-2*m_plot-(ny-1)*s_plot) / ny;
    x_pos = repmat(m_plot+(0:(nx-1))*w_plot+(0:(nx-1))*s_plot, ny, 1)';
    y_pos = repmat(1-(h_plot+m_plot+(0:(ny-1))*h_plot+(0:(ny-1))*s_plot), nx, 1);
    idxOut = ones(size(WP.river));
    fig1 = figure(); 
    subplot('Position', [x_pos(1) y_pos(1) h_plot w_plot]); hold on;
        plotWPdata(WP.So, WP.Cm_ ./ WP.So, categorical(WP.river), idxOut, fig1);
        plot(fivepctexceed.slope, fivepctexceed.cms0, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'k')
        text(0.9, 0.9, 'a', 'units', 'normalized', 'FontSize', 14)
        annotation('arrow', [0.385 0.385], [0.625, 0.755])
        text(0.72, 0.64, {'increasing','discharge'}, 'units', 'normalized', 'FontSize', 10)
        set(gca, 'XScale', 'log')
        xlabel('$S_0$')
        ylabel('$\bar{C} / S_0$')   
        xlim([1e-5 1e-2])
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot('Position', [x_pos(2) y_pos(2) h_plot w_plot]); hold on;
        plotWPdata(WP.So, WP.vs50_u_, categorical(WP.river), idxOut, fig1);
        plot(fivepctexceed.slope, fivepctexceed.uws, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'k')
        text(0.9, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
        annotation('arrow', [0.865 0.865], [0.755, 0.625])
        text(0.72, 0.1, {'increasing','discharge'}, 'units', 'normalized', 'FontSize', 10)
        set(gca, 'XScale', 'log')
        xlabel('$S_0$')
        ylabel('$w_{s,50}/u_*$')
        xlim([1e-5 1e-2])
        ylim([0 1])
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot('Position', [x_pos(3) y_pos(3) h_plot w_plot]); hold on;
        [WPmarkers] = plotWPdata(WP.So, WP.Kred, categorical(WP.river), idxOut, fig1);
        plot(fivepctexceed.slope, fivepctexceed.kred, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'k')
        text(0.9, 0.9, 'c', 'units', 'normalized', 'FontSize', 14)
        annotation('arrow', [0.385 0.385], [0.28, 0.150])
        text(0.72, 0.1, {'increasing','discharge'}, 'units', 'normalized', 'FontSize', 10)
        set(gca, 'XScale', 'log')
        xlabel('$S_0$')
        ylabel('$\overline{K}_{red}$')
        xlim([1e-5 1e-2])
        ylim([0 1])
        axis square
        box on
        lgd = legend(WPmarkers, cellstr(unique(WP.river)), 'location', 'southwest', 'Interpreter', 'latex');
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    subplot('Position', [x_pos(4) y_pos(4) h_plot w_plot]); hold on;
        plotWPdata(WP.C5t ./ WP.So, WP.Kred, categorical(WP.river), idxOut, fig1);
        plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'kx-', 'LineWidth', 1)
        text(0.9, 0.9, 'd', 'units', 'normalized', 'FontSize', 14)
        ylabel('$\alpha_{WP04} \equiv \overline{K}_{red}$')
        xlabel('$\bar{c}_b / S_0$')
        xlim([0 100])
        ylim([0 1])
        axis square
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    set(gcf, 'Pos', [100 100 600 600], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    if printOut
        print('-dpng', '-r300', [print_root, 'WP_plots_remake.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'WP_plots_remake.eps']);
        pause(1);
    end

    %% every concProf
    if true
        [sortByUstar, sortByUstarIdx] = sort( cell2mat(arrayfun(@(s) s.Velocity.ustarCalib, df, 'Unif', 0)) ); %#ok<*ASGLU>
        plotMatLen = ceil(sqrt(size(df,1)));
        m_plot = 0.025; % margins
        s_plot = 0.0125; % spacing
        w_plot = (1-2*m_plot-6*s_plot) / 7;
        h_plot = (1-2*m_plot-7*s_plot) / 8;
        x_pos = repmat(m_plot+(0:6)*w_plot+(0:6)*s_plot, 8, 1)';
        y_pos = repmat(1-(h_plot+m_plot+(0:7)*h_plot+(0:7)*s_plot), 7, 1);
        figure()
        for p = 1:size(df,1)
            dfIdx = sortByUstarIdx(p);
            if datenum(df(dfIdx).CollectionDate) < datenum('01/01/2016') 
                yrCol = col1;
                yrMark = mar1;
            elseif datenum(df(dfIdx).CollectionDate) < datenum('01/01/2017')
                yrCol = col2;
                yrMark = mar2;
            else
                yrCol = col3;
                yrMark = mar3;
            end           
            subplot('Position', [x_pos(p) y_pos(p) w_plot, h_plot]); hold on;
                [gsclassfit] = plot(df(dfIdx).concProf.gsClassModel.CsSum, df(dfIdx).concProf.gsClassModel.ZsSum, '-', 'Color', [0.6 0.6 0.6], 'MarkerSize', 3, 'LineWidth', 1.5);
                [gsclasspred] = plot(df(dfIdx).concProf.gsClassPred.CsSum, df(dfIdx).concProf.gsClassPred.ZsSum, 'k-', 'MarkerSize', 3, 'LineWidth', 1.5);
                [wpalpha] = plot(df(dfIdx).concProf.WPalphaPred.CsSum, df(dfIdx).concProf.WPalphaPred.ZsSum, ':', 'Color', [0.482 0.3059 0.5843], 'MarkerSize', 3, 'LineWidth', 2);
                [MYsum_mdl] = plot(df(dfIdx).concProf.MYfullPred.CsSum, df(dfIdx).concProf.MYfullPred.ZsSum, ...
                    '--', 'Color', [0.369 0.71 0.349], 'LineWidth', 1.5); % mellor yamada prediction
                [samps] = plot(df(dfIdx).waterSamplesTable.concNW, df(dfIdx).waterSamplesTable.sampleZ, 'k', ...
                    'MarkerFaceColor', yrCol, 'MarkerSize', 4, 'Marker', yrMark, 'LineStyle', 'none');
                box on; set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', ...
                    'layer', 'top', 'TickLabelInterpreter', 'latex')
                xl = [0 max(df(dfIdx).waterSamplesTable.concNW)*1.3];
                xlim(xl);
                yl = [0 df(dfIdx).FlowDepthAtCollection*1.2];
                ylim(yl);
                set(gca, 'xTick', (0:floor(xl(end))))
                set(gca, 'yTick', (0:0.5:floor(yl(end))))
                xText = 0.45;
                text(xText, 0.85, strrep(df(dfIdx).StationID, '_', ' '), 'FontSize', 9, 'Units', 'normalized')
                text(xText, 0.75, ['$u_*$ = ', num2str(round(df(dfIdx).Velocity.ustarCalib*100,2)), ' cm/s'], 'FontSize', 8, 'Units', 'normalized')
                text(xText, 0.65, ['$c_b$ = ', num2str(round(df(dfIdx).concProf.totalModel.Cs(1),2)), ' g/L'], 'FontSize', 8, 'Units', 'normalized')
                text(xText, 0.55, ['H = ', num2str(round(df(dfIdx).FlowDepthAtCollection,1)), ' m'], 'FontSize', 8, 'Units', 'normalized')
                set(gca, 'xticklabels', [])
                set(gca, 'yticklabels', [])
                set(gca,'TickLength',[0.015, 0.01])
        end
        lg = legend([gsclassfit gsclasspred wpalpha MYsum_mdl samps], {'$\alpha_f$ profile', '$\alpha_{1.0}$ profile', ...
            '$\alpha_{WP04} profile$', '$MY$ profile', ...
            'no-washload samples'}, 'Interpreter', 'latex', 'EdgeColor', [0 0 0]);
        p = p + 1;
        set(lg, 'Position', [x_pos(p), y_pos(p)+s_plot, w_plot-s_plot, h_plot-4*s_plot])
        p = p + 1;
        subplot('Position', [x_pos(p) y_pos(p) w_plot, h_plot])
            xlabel('concentration (g/L)')
            ylabel('dist. above bed (m)')
            box on; set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gca, 'xticklabels', [])
            set(gca, 'yticklabels', [])
            set(gca,'TickLength',[0.015, 0.01])
        set(gcf, 'Pos', [50 100 1600 1500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        if printOut
            print('-dpng', '-r200', [print_root, 'everyConcProf.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'everyConcProf.eps']);
            pause(1);
        end
    end
    
    if true
        velProfIdx = arrayfun(@(x) ~(size(fieldnames(x.velProf),1)<5), df);
        velProfIdxLocs = find(velProfIdx);
        figure();
        for p = 1:sum(velProfIdx)
            subplot(3, 5, p); hold on;
            dfIdx = velProfIdxLocs(p);
            if datenum(df(dfIdx).CollectionDate) < datenum('01/01/2016') 
                yrCol = col1;
                yrMark = mar1;
            elseif datenum(df(dfIdx).CollectionDate) < datenum('01/01/2017')
                yrCol = col2;
                yrMark = mar2;
            else
                yrCol = col3;
                yrMark = mar3;
            end
            if ~isempty(df(dfIdx).Velocity.adcp.mean)
                errorbar(df(dfIdx).Velocity.adcpTable.mean, df(dfIdx).Velocity.adcpTable.measZ, df(dfIdx).Velocity.adcpTable.std, 'k.', 'horizontal')
                plot(df(dfIdx).Velocity.adcpTable.mean, df(dfIdx).Velocity.adcpTable.measZ, 'k', 'MarkerFaceColor', yrCol, 'MarkerSize', 5, 'Marker', yrMark, 'LineStyle', 'none')
                plot(df(dfIdx).velProf.modelNoAlpha.Us, df(dfIdx).velProf.modelNoAlpha.Zs, 'k-')
            end
            if ~isempty(df(dfIdx).Velocity.meter.velmean)
                if ~iscell(df(dfIdx).Velocity.meterTable.err(1))
                    errorbar(df(dfIdx).Velocity.meterTable.mean, df(dfIdx).Velocity.meterTable.measZ, df(dfIdx).Velocity.meterTable.err, 'k.', 'horizontal')
                else
                    df(dfIdx).Velocity.meterTable.errneg = df(dfIdx).Velocity.meterTable.mean - arrayfun(@(x) x{1}(1), df(dfIdx).Velocity.meterTable.err, 'Unif', 1);
                    df(dfIdx).Velocity.meterTable.errpos = arrayfun(@(x) x{1}(2), df(dfIdx).Velocity.meterTable.err, 'Unif', 1) - df(dfIdx).Velocity.meterTable.mean;
                    errorbar(df(dfIdx).Velocity.meterTable.mean, df(dfIdx).Velocity.meterTable.measZ, df(dfIdx).Velocity.meterTable.errneg, df(dfIdx).Velocity.meterTable.errpos, 'k.', 'horizontal')
                end
                plot(df(dfIdx).Velocity.meterTable.mean, df(dfIdx).Velocity.meterTable.measZ, 'k', 'MarkerSize', 5, 'Marker', yrMark, 'LineStyle', 'none')
                plot(df(dfIdx).velProf.modelKred.Us, df(dfIdx).velProf.modelKred.Zs, 'k-')
            end
            params1a = {strcat(strrep(df(dfIdx).StationID, '_', '\_')), ...
                strcat('$z_0 = ', sprintsci(df(dfIdx).velProf.modelNoAlpha.params.z0), '$'), ...
                strcat('$u_* = ', num2str(round(df(dfIdx).velProf.modelNoAlpha.params.ustar,4)), '$'), ...
                strcat('$R^2 = ', num2str(round(df(dfIdx).velProf.modelNoAlpha.params.model.Rsquared.Ordinary,2)), '$')};
            text(0.05, 0.62, sprintf('%s\n', params1a{:}), 'Parent', gca, 'units', 'normalized');
            ylim([0 df(dfIdx).FlowDepthAtCollection*1.1])
            xlim([0 max(df(dfIdx).velProf.modelNoAlpha.Us)*1.5])
            xlabel('velocity (m/s)')
            ylabel('distance above bed (m)')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 1100 800], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        end
        if printOut
            print('-dpng', '-r300', [print_root, 'everyVelProf.png']);
            print('-depsc', '-painters', '-r300', [print_root, '/everyVelProf.eps']);
            pause(1);
        end
    
    end

    %% examples and summary figures
    m_plot = 0.075; % margins
    s_plot = 0.075; % spacing
    nx = 3;
    ny = 3;
    w_plot = (1-2*m_plot-(nx-1)*s_plot) / nx;
    h_plot = (1-2*m_plot+0.5*m_plot-(ny-1)*s_plot) / ny;
    x_pos = repmat(m_plot+(0:(nx-1))*w_plot+(0:(nx-1))*s_plot, ny, 1)';
    y_pos = repmat(1-(h_plot+1.5*m_plot+(0:(ny-1))*h_plot+(0:(ny-1))*s_plot), nx, 1);
    exNums = [27 15 3]; % df idx to plot for 2015, 2016, 2018 data
    figure()
    for p = 1:length(exNums)
        if datenum(df(exNums(p)).CollectionDate) < datenum('01/01/2016') 
            yrCol = col1;
            yrMark = mar1;
        elseif datenum(df(exNums(p)).CollectionDate) < datenum('01/01/2017')
            yrCol = col2;
            yrMark = mar2;
        else
            yrCol = col3;
            yrMark = mar3;
        end      
        subplot('position', [x_pos(p) y_pos(p) w_plot, h_plot]); hold on;
            % gs class specific models in background
            plot(df(exNums(p)).concProf.gsClassModel.Cs, df(exNums(p)).concProf.gsClassModel.Zs, '--', 'Color', [0.6 0.6 0.6])
            plot(concatenate_cell_vectors(df(exNums(p)).waterSamplesTable.concNWbyClass), repelem(df(exNums(p)).waterSamplesTable.sampleZ, length(df(exNums(p)).waterSamplesTable.gsClass{1})), ...
                'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'Marker', yrMark, 'LineStyle', 'none')
            gsclassModel_lines = plot(NaN, NaN, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', ...
                'MarkerSize', 3, 'Marker', yrMark, 'LineStyle', '--', 'Color', [0.6 0.6 0.6]); % blank version, for the legend
            
            % main samples and models
            sumGsClass_fit = plot(df(exNums(p)).concProf.gsClassModel.CsSum, df(exNums(p)).concProf.gsClassModel.ZsSum, ...
                '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5); % summed class fits
            sumGsClass_mdl = plot(df(exNums(p)).concProf.gsClassPred.CsSum, df(exNums(p)).concProf.gsClassPred.ZsSum, ...
                'k-', 'LineWidth', 1.5); % sum of predictions
            [wpalpha] = plot(df(exNums(p)).concProf.WPalphaPred.CsSum, df(exNums(p)).concProf.WPalphaPred.ZsSum, ':', ...
               'Color', [0.482 0.3059 0.5843], 'MarkerSize', 3, 'LineWidth', 2); % prediction from WPalpha
            MYsum_mdl = plot(df(exNums(p)).concProf.MYfullPred.CsSum, df(exNums(p)).concProf.MYfullPred.ZsSum, ...
                '--', 'Color', [0.369 0.71 0.349], 'LineWidth', 1.5); % wright parker prediction
            NW_samples = plot(df(exNums(p)).waterSamplesTable.concNW, df(exNums(p)).waterSamplesTable.sampleZ, ...
                'k', 'MarkerFaceColor', yrCol, 'MarkerSize', 6, 'Marker', yrMark, 'LineStyle', 'none'); % NW samples
            
            
            % full measured samples
            wW_samples = plot(df(exNums(p)).waterSamplesTable.conc, df(exNums(p)).waterSamplesTable.sampleZ, ...
                'k', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 4, 'Marker', yrMark, 'LineStyle', 'none');
            
            box on; set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            xlabel('conc (g/L)')
            ylabel('dist. above bed (m)')
            ylim([0 df(exNums(p)).FlowDepthAtCollection])
            text(0.9, 0.9, char('a'+p-1), 'units', 'normalized', 'FontSize', 14)
    end
    lg = legend([NW_samples, wW_samples, gsclassModel_lines, sumGsClass_fit, sumGsClass_mdl, MYsum_mdl, wpalpha], ...
            {'no-washload samples', 'raw samples', '$\alpha_{f,i}$ profiles', '$\alpha_f$ profile', ...
            '$\alpha_{1.0}$ profile', '$MY$ profile', '$\alpha_{WP04}$ profile'}, 'Interpreter', 'latex');
    set(lg, 'Position', [0.1 0.9 0.8 0.1])
    lg.NumColumns = 3;
    bw = 0.05;
    subplot('position', [x_pos(4) y_pos(4) w_plot, h_plot]); hold on;
        plot(fliplr(jt.waterTable15ConcNWWide), jt.waterTable15WideJitter, 'LineStyle', 'none', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
        boxplot(fliplr(jt.waterTable15ConcNWWide), 'positions', waterSampleZs, 'orientation', 'horizontal', 'widths', bw)
        text(0.9, 0.9, 'd', 'units', 'normalized', 'FontSize', 14)
            xlim([0 50])
            ylim([-0.1 1])
            ylabel('height above bed $(z/H)$')
            xlabel('concentration (g/L)')
            set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.50','0.9'})
            cleanup_boxplot()
            pause(0.2)
    subplot('position', [x_pos(5) y_pos(5) w_plot, h_plot]); hold on;
        plot(fliplr(jt.waterTable16ConcNWWide), jt.waterTable16WideJitter, 'LineStyle', 'none', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
        boxplot(fliplr(jt.waterTable16ConcNWWide), 'positions', waterSampleZs, 'orientation', 'horizontal', 'widths', bw)
        text(0.9, 0.9, 'e', 'units', 'normalized', 'FontSize', 14)
            xlim([0 50])
            ylim([-0.1 1])
            ylabel('height above bed $(z/H)$')
            xlabel('concentration (g/L)')
            set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.50','0.9'})
            cleanup_boxplot()
            pause(0.2)
    subplot('position', [x_pos(6) y_pos(6) w_plot, h_plot]); hold on;
        plot(fliplr(jt.waterTable18ConcNWWide), jt.waterTable18WideJitter, 'LineStyle', 'none', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
        boxplot(fliplr(jt.waterTable18ConcNWWide), 'positions', waterSampleZs, 'orientation', 'horizontal', 'widths', bw)
        text(0.9, 0.9, 'f', 'units', 'normalized', 'FontSize', 14)
            xlim([0 50])
            ylim([-0.1 1])
            ylabel('height above bed $(z/H)$')
            xlabel('concentration (g/L)')
            set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.50','0.9'})
            cleanup_boxplot()
            pause(0.2)
    subplot('position', [x_pos(7) y_pos(7) w_plot, h_plot]); hold on;
        plot(fliplr(jt.samplesTable15d5NWWide), jt.samplesTable15WideJitter, 'LineStyle', 'none', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
        plot(fliplr(jt.samplesTable15d50NWWide),jt.samplesTable15WideJitter, 'LineStyle', 'none', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
        plot(fliplr(jt.samplesTable15d90NWWide), jt.samplesTable15WideJitter, 'LineStyle', 'none', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
        boxplot(fliplr(jt.samplesTable15d5NWWide), 'positions', sampleZs, 'orientation', 'horizontal', 'widths', bw)
        boxplot(fliplr(jt.samplesTable15d50NWWide), 'positions', sampleZs, 'orientation', 'horizontal', 'widths', bw)
        boxplot(fliplr(jt.samplesTable15d90NWWide), 'positions', sampleZs, 'orientation', 'horizontal', 'widths', bw)
        text(0.9, 0.9, 'g', 'units', 'normalized', 'FontSize', 14)
            xlim([0 200])
            ylim([-0.1 1])
            ylabel('height above bed $(z/H)$')
            xlabel('$D_{5}$, $D_{50}$, $D_{90}$ grain size $(\mu$m$)$')
            set(gca, 'YTick', sampleZs, 'YTickLabel', {'bed', '0.05','0.15','0.25','0.50','0.9'})
            cleanup_boxplot()
            pause(0.2)
    subplot('position', [x_pos(8) y_pos(8) w_plot, h_plot]); hold on;
        plot(fliplr(jt.samplesTable16d5NWWide), jt.samplesTable16WideJitter, 'LineStyle', 'none', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
        plot(fliplr(jt.samplesTable16d50NWWide),jt.samplesTable16WideJitter, 'LineStyle', 'none', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
        plot(fliplr(jt.samplesTable16d90NWWide), jt.samplesTable16WideJitter, 'LineStyle', 'none', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
        boxplot(fliplr(jt.samplesTable16d5NWWide), 'positions', sampleZs, 'orientation', 'horizontal', 'widths', bw)
        boxplot(fliplr(jt.samplesTable16d50NWWide), 'positions', sampleZs, 'orientation', 'horizontal', 'widths', bw)
        boxplot(fliplr(jt.samplesTable16d90NWWide), 'positions', sampleZs, 'orientation', 'horizontal', 'widths', bw)
        text(0.9, 0.9, 'h', 'units', 'normalized', 'FontSize', 14)
            xlim([0 200])
            ylim([-0.1 1])
            ylabel('height above bed $(z/H)$')
            xlabel('$D_{5}$, $D_{50}$, $D_{90}$ grain size $(\mu$m$)$')
            set(gca, 'YTick', sampleZs, 'YTickLabel', {'bed', '0.05','0.15','0.25','0.50','0.9'})
            cleanup_boxplot()
            pause(0.2)
    subplot('position', [x_pos(9) y_pos(9) w_plot, h_plot]); hold on;
        plot(fliplr(jt.samplesTable18d5NWWide), jt.samplesTable18WideJitter, 'LineStyle', 'none', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
        plot(fliplr(jt.samplesTable18d50NWWide), jt.samplesTable18WideJitter, 'LineStyle', 'none', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
        plot(fliplr(jt.samplesTable18d90NWWide), jt.samplesTable18WideJitter, 'LineStyle', 'none', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
        boxplot(fliplr(jt.samplesTable18d5NWWide), 'positions', sampleZs, 'orientation', 'horizontal', 'widths', bw)
        boxplot(fliplr(jt.samplesTable18d50NWWide), 'positions', sampleZs, 'orientation', 'horizontal', 'widths', bw)
        boxplot(fliplr(jt.samplesTable18d90NWWide), 'positions', sampleZs, 'orientation', 'horizontal', 'widths', bw)
        text(0.9, 0.9, 'i', 'units', 'normalized', 'FontSize', 14)
            xlim([0 200])
            ylim([-0.1 1])
            ylabel('height above bed $(z/H)$')
            xlabel('$D_{5}$, $D_{50}$, $D_{90}$ grain size $(\mu$m$)$')
            set(gca, 'YTick', sampleZs, 'YTickLabel', {'bed', '0.05','0.15','0.25','0.50','0.9'})
            cleanup_boxplot()
            pause(0.2)
    set(gcf, 'Pos', [100 100 700 700]);
    set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    drawnow; pause(0.7);
    if printOut
        print('-dpng', '-r300', [print_root, 'threeYearsExample.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'threeYearsExample.eps']);
        pause(1);
    end
    
    
    %% washload percentage for supplement
    figure();
        subplot(1, 3, 1); hold on;
            plot(fliplr(jt.waterTable15WLpercentWide), jt.waterTable15WideJitter, 'o', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
            boxplot(fliplr(jt.waterTable15WLpercentWide), 'positions', waterSampleZs, 'orientation', 'horizontal')
                xlim([0 100])
                ylabel('height above bed $(z/H)$')
                xlabel('washload percentage (\%)')
                set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.5','0.9'})
                cleanup_boxplot()
        subplot(1, 3, 2); hold on;
            plot(fliplr(jt.waterTable16WLpercentWide), jt.waterTable16WideJitter, 'o', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
            boxplot(fliplr(jt.waterTable16WLpercentWide), 'positions', waterSampleZs, 'orientation', 'horizontal')
                xlim([0 100])
                ylabel('height above bed $(z/H)$')
                xlabel('washload percentage (\%)')
                set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.5','0.9'})
                cleanup_boxplot()
        subplot(1, 3, 3); hold on;
            plot(fliplr(jt.waterTable18WLpercentWide), jt.waterTable18WideJitter, 'o', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
            boxplot(fliplr(jt.waterTable18WLpercentWide), 'positions', waterSampleZs, 'orientation', 'horizontal')
                xlim([0 100])
                ylabel('height above bed $(z/H)$')
                xlabel('washload percentage (\%)')
                set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.5','0.9'})
                cleanup_boxplot()
        set(gcf, 'Pos', [100 100 1050 350]);
        set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        drawnow; pause(0.1);
        if printOut
            print('-dpng', '-r300', [print_root, 'washloadfraction.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'washloadfraction.eps']);
            pause(1);
        end
    
    
    %% grainsize overbank for supplement
    [D10PredJitterValues, DjitterLocs, YrMat] = categorical2jitterMat(station.transport.pred_D10, arrayfun(@(x) year(x.CollectionDate), df));
    [D50PredJitterValues, ~] = categorical2jitterMat(station.transport.pred_D50, arrayfun(@(x) year(x.CollectionDate), df));
    [D90PredJitterValues, ~] = categorical2jitterMat(station.transport.pred_D90, arrayfun(@(x) year(x.CollectionDate), df));
    [D10FitJitterValues, ~] = categorical2jitterMat(station.transport.fit_D10, arrayfun(@(x) year(x.CollectionDate), df));
    [D50FitJitterValues, ~] = categorical2jitterMat(station.transport.fit_D50, arrayfun(@(x) year(x.CollectionDate), df));
    [D90FitJitterValues, ~] = categorical2jitterMat(station.transport.fit_D90, arrayfun(@(x) year(x.CollectionDate), df));
    YrMat = YrMat(:);

    DjitterLocs = [DjitterLocs(:,1)-1; DjitterLocs(:,2)-2; DjitterLocs(:,3)-3];
    
    f = figure();
        plotYRdata(DjitterLocs + 1, D10PredJitterValues(:), YrMat, ones(size(YrMat,1),1), f, 'ColorEdgeOnly', true)
        plotYRdata(DjitterLocs + 2, D10PredJitterValues(:), YrMat, ones(size(YrMat,1),1), f, 'ColorEdgeOnly', true)
        plotYRdata(DjitterLocs + 4, D50PredJitterValues(:), YrMat, ones(size(YrMat,1),1), f, 'ColorEdgeOnly', true)
        plotYRdata(DjitterLocs + 5, D50PredJitterValues(:), YrMat, ones(size(YrMat,1),1), f, 'ColorEdgeOnly', true)
        plotYRdata(DjitterLocs + 7, D90PredJitterValues(:), YrMat, ones(size(YrMat,1),1), f, 'ColorEdgeOnly', true)
        plotYRdata(DjitterLocs + 8, D90PredJitterValues(:), YrMat, ones(size(YrMat,1),1), f, 'ColorEdgeOnly', true)
        boxplot([station.transport.pred_D10, station.transport.fit_D10,...
             station.transport.pred_D50, station.transport.fit_D50,...
             station.transport.pred_D90, station.transport.fit_D90], 'positions', [1 2 4 5 7 8]); 
        set(gca, 'Xtick', [1 2 4 5 7 8], 'XTickLabels', {'$D_{10,1.0}$', '$D_{10,f}$', '$D_{50,1.0}$', ...
                                                         '$D_{50,f}$', '$D_{90,1.0}$', '$D_{90,f}$'});
        ylabel('grain size ($\mu m$)')
        cleanup_boxplot()

    set(gcf, 'Pos', [100 100 500 350]);
    set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    drawnow; pause(0.1);
    if printOut
        print('-dpng', '-r300', [print_root, 'grainsize_overbank.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'grainsize_overbank.eps']);
        pause(1);
    end

    
    %% single large example plot for talks
    figure(); hold on;
    if datenum(df(exNums(p)).CollectionDate) < datenum('01/01/2016') 
        yrCol = col1;
        yrMark = mar1;
    elseif datenum(df(exNums(p)).CollectionDate) < datenum('01/01/2017')
        yrCol = col2;
        yrMark = mar2;
    else
        yrCol = col3;
        yrMark = mar3;
    end      
    % gs class specific models in background
    plot(df(exNums(p)).concProf.gsClassModel.Cs, df(exNums(p)).concProf.gsClassModel.Zs, '--', 'Color', [0.6 0.6 0.6])
    plot(concatenate_cell_vectors(df(exNums(p)).waterSamplesTable.concNWbyClass), repelem(df(exNums(p)).waterSamplesTable.sampleZ, length(df(exNums(p)).waterSamplesTable.gsClass{1})), ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'Marker', yrMark, 'LineStyle', 'none')
    gsclassModel_lines = plot(NaN, NaN, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 3, 'Marker', yrMark, 'LineStyle', '--', 'Color', [0.6 0.6 0.6]); % blank version, for the legend

    % main samples and models
    sumGsClass_fit = plot(df(exNums(p)).concProf.gsClassModel.CsSum, df(exNums(p)).concProf.gsClassModel.ZsSum, ...
        '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5); % summed class fits
    sumGsClass_mdl = plot(df(exNums(p)).concProf.gsClassPred.CsSum, df(exNums(p)).concProf.gsClassPred.ZsSum, ...
        'k-', 'LineWidth', 1.5); % sum of predictions
    [wpalpha] = plot([NaN], [NaN], 'color', 'none');
    MYsum_mdl = plot(df(exNums(p)).concProf.MYfullPred.CsSum, df(exNums(p)).concProf.MYfullPred.ZsSum, ...
        '--', 'Color', [0.369 0.71 0.349], 'LineWidth', 1.5); % wright parker prediction
    NW_samples = plot(df(exNums(p)).waterSamplesTable.concNW, df(exNums(p)).waterSamplesTable.sampleZ, ...
        'k', 'MarkerFaceColor', yrCol, 'MarkerSize', 6, 'Marker', yrMark, 'LineStyle', 'none'); % NW samples


    % full measured samples
    wW_samples = plot(df(exNums(p)).waterSamplesTable.conc, df(exNums(p)).waterSamplesTable.sampleZ, ...
        'k', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 4, 'Marker', yrMark, 'LineStyle', 'none');

    box on; set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    xlabel('conc (g/L)')
    ylabel('dist. above bed (m)')
    ylim([0 df(exNums(p)).FlowDepthAtCollection])
    text(0.9, 0.9, char('a'+p-1), 'units', 'normalized', 'FontSize', 14)
    lg = legend([wW_samples, NW_samples, gsclassModel_lines, sumGsClass_fit, sumGsClass_mdl, MYsum_mdl, wpalpha], ...
            {'raw samples', 'no-washload samples', '$\alpha_{f,i}$ profiles', '$\alpha_f$ profile', ...
            '$\alpha_{1.0}$ profile', '$MY$ profile', '$\alpha_{WP04}$ profile'}, 'Interpreter', 'latex');
    set(gcf, 'Pos', [100 100 400 400]);
    set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    drawnow; pause(0.7);
    if printOut
        print('-dpng', '-r300', [print_root, 'largeSingleExample.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'largeSingleExample.eps']);
        pause(1);
    end
    
    
    %% make the supplementary grid random search plot
    load('stratificationParameterSpaceRand.mat', 'out', 'K_red', 'D_rand', 'cb_rand', 'S_rand', 'h_rand')
    ustar_rand = sqrt(9.81.*h_rand.*S_rand);
    rs = array2table([cb_rand, S_rand, D_rand, h_rand, ustar_rand, K_red(:)], ...
                     'VariableNames', {'cb', 'S', 'D', 'h', 'ustar', 'K_red'}); 
    rs.ws = get_DSV(rs.D, 0.7, 3.5, con);
    rs.ustarws = rs.ustar ./ rs.ws;
    rs.Rep = (sqrt(con.g .* con.R .* (rs.D)) .* rs.D) ./ con.nu;
    rs.Dstar = ((con.R * con.g) / (con.nu * con.nu))^(1/3) .* rs.D;
    rs.logS = log10(rs.S);
    rs.logcb = log10(rs.cb);
    f= figure(); 
    subplot(1, 2, 1); hold on;
        scatter(rs.cb ./ rs.S, rs.K_red, 55, rs.D, 'filled', 'MarkerFaceAlpha', 0.5);
        fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 4, 'LineWidth', 2)
        set(gca, 'XScale', 'log')
        ylim([0, 1])
        xlim([1e-1, 1e2])
        xlabel('$c_b / S$')
        ylabel('$\alpha \equiv \overline{K}_{red}$')
        colb = colorbar;
        title(colb, '$D$', 'interpreter', 'latex');
        colb.TickLabelInterpreter = 'latex';
        caxis([0 500e-6])
        box on
        legend(WP_marks, unique(WP.river), 'location', 'southwest', 'Interpreter', 'latex')
        box on; set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    xx = logspace(-1, 2, 50);
    yy = linspace(0, 1, 50);
    subplot(1, 2, 2); hold on;
        h = histogram2(rs.cb ./ rs.S, rs.K_red, xx,yy, 'FaceColor','flat', 'EdgeColor', 'none');
        plot3(mods.WPalphaEvalXs, mods.WPalphaEval, 100*ones(size(mods.WPalphaEval)), 'k-x', 'MarkerSize', 4, 'LineWidth', 2)
        set(gca, 'XScale', 'log')
        colb = colorbar;
        title(colb, 'count density', 'interpreter', 'latex');
        colb.TickLabelInterpreter = 'latex';
        colormap(gca, flipud(hot))
        set(gca, 'colorscale', 'linear')
        ylim([0, 1])
        xlim([1e-1, 1e2])
        xlabel('$c_b / S$')
        ylabel('$\alpha \equiv \overline{K}_{red}$')
        box on; set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
    set(gcf, 'Pos', [100 100 800 350]);
    set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    drawnow; pause(0.5);
    if printOut
        print('-dpng', '-r300', [print_root, 'randomparamspace.png']);
        print('-depsc', '-painters', '-r300', [print_root, 'randomparamspace.eps']);
        pause(1);
    end
    
    %% here are all the variables that get written into the text, to update if something changes
    
    % do grain size changes!
    
    % entrainment data trends
    entrtab = array2table([mean(eAlpha1(low_idx)), mean(eAlpha1(in_idx)), mean(eAlpha1(high_idx));...
                           std(eAlpha1(low_idx)),  std(eAlpha1(in_idx)),  std(eAlpha1(high_idx));...
                           median(eAlpha1(low_idx)),  median(eAlpha1(in_idx)),  median(eAlpha1(high_idx));...
                           NaN, ranksum(eAlpha1(low_idx), eAlpha1(in_idx)),  ranksum(eAlpha1(low_idx), eAlpha1(high_idx));...
                           ranksum(eAlpha1(in_idx), eAlpha1(low_idx)), NaN, ranksum(eAlpha1(in_idx), eAlpha1(high_idx));...
                           ranksum(eAlpha1(high_idx), eAlpha1(low_idx)),  ranksum(eAlpha1(high_idx), eAlpha1(in_idx)), NaN], 'VariableNames', {'low', 'in', 'high'});
    
    
    % slope regression
    mfunct_str = [num2str(rndw), ' * log10(S) + ', num2str(rndx)];
    nfunct_str = [num2str(rndy), ' * log10(S) + ', num2str(rndz)];
    slope_mdl_r2 = slope_mdl_r2; 
    
    % grain size i alpha regression
    regmdl_str = [num2str(arnd), '[u_*/ws cb5/S]^(', num2str(brnd), ')'];
    reg_mdl_r2 = reg_mdl_r2; 
    
    
    % overbank sediment dicussion section
    obtests  = [signrank(station.transport.pred_D10,station.transport.fit_D10),...
                signrank(station.transport.pred_D50,station.transport.fit_D50),...
                signrank(station.transport.pred_D90,station.transport.fit_D90)];
    obdata   = [nanmean((station.transport.pred_D10-station.transport.fit_D10)./station.transport.pred_D10),...
                nanmean((station.transport.pred_D50-station.transport.fit_D50)./station.transport.pred_D50),...
                nanmean((station.transport.pred_D90-station.transport.fit_D90)./station.transport.pred_D90)];
    
end

function [SSE] = minimize_regress(x, expected_vals, varargin)
    % x is vector containing a, b, c, d; coeffs on m = f(S) and n = f(S)
    % varargin is dependent on the prediction function: u*/wsRep^0.9, and S
    
    a = x(1);
    b = x(2);
    c = x(3);
    d = x(4);
    
    m = a .* log10(varargin{2}) + b;
    n = c .* log10(varargin{2}) + d;
    
    predicted_vals = 1 - (m .* varargin{1} .^ n);
    
    y = expected_vals;
    yhat = predicted_vals;
    r = y - yhat;
    SSE = norm(r, 2).^2;

end

function [SSE] = minimize_regress2(x, expected_vals, varargin)
    % x is vector containing a, b, c; coeffs on m = f(S) and n = f(S)
    % varargin is dependent on the prediction function: u*/wsRep^0.?, and S
    
    a = x(1);
    b = x(2);
    c = x(3);
    m = a .* log10(varargin{2}) + b;
    n = c;
    
    predicted_vals = 1 - (m .* varargin{1} .^ n);
    
    y = expected_vals;
    yhat = predicted_vals;
    r = y - yhat;
    SSE = norm(r, 2).^2;

end
