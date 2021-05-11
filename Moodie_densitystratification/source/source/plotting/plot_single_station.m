function plot_single_station(df)
    
    global col1 col2 col3

    printOut = true;
    
    figure()
    
        %% velocity profiles
        subplot(2, 2, 1); hold on
            % handle adcp if exists
            if ~isempty(df.Velocity.adcp)
                errorbar(df.Velocity.adcpTable.mean, df.Velocity.adcpTable.measZ, df.Velocity.adcpTable.std, 'k.', 'horizontal')
                plot(df.Velocity.adcpTable.mean, df.Velocity.adcpTable.measZ, 'ko', 'MarkerFaceColor', col3, 'MarkerSize', 5)
                plot(df.velProf.modelAdcpNoAlpha.Us, df.velProf.modelAdcpNoAlpha.Zs, 'k-', 'LineWidth', 1.5)
                params1b = {strcat('z_0 = ', sprintsci(df.velProf.modelAdcpNoAlpha.params.z0)), ...
                strcat('u_* = ', num2str(round(df.velProf.modelAdcpNoAlpha.params.ustar,4))), ...
                strcat('R^2 = ', num2str(round(df.velProf.modelAdcpNoAlpha.params.model.Rsquared.Ordinary,2)))};
                text(0.6, 0.7, sprintf('%s\n', params1b{:}), 'Parent', gca, 'units', 'normalized');
            end
            % handle meter if exists
            if ~isempty(df.Velocity.meter)
                if ~iscell(df.Velocity.meterTable.err(1))
                    errorbar(df.Velocity.meterTable.mean, df.Velocity.meterTable.measZ, df.Velocity.meterTable.err, 'k.', 'horizontal')
                else
                    df.Velocity.meterTable.errneg = df.Velocity.meterTable.mean - arrayfun(@(x) x{1}(1), df.Velocity.meterTable.err, 'Unif', 1);
                    df.Velocity.meterTable.errpos = arrayfun(@(x) x{1}(2), df.Velocity.meterTable.err, 'Unif', 1) - df.Velocity.meterTable.mean;
                    errorbar(df.Velocity.meterTable.mean, df.Velocity.meterTable.measZ, df.Velocity.meterTable.errneg, df.Velocity.meterTable.errpos, 'k.', 'horizontal')
                end
                plot(df.Velocity.meterTable.mean, df.Velocity.meterTable.measZ, 'k^', 'MarkerFaceColor', col3, 'MarkerSize', 5)
                plot(df.velProf.modelMeterNoAlpha.Us, df.velProf.modelMeterNoAlpha.Zs, 'k-', 'LineWidth', 1.5)
                params1a = {strcat('z_0 = ', sprintsci(df.velProf.modelMeterNoAlpha.params.z0)), ...
                strcat('u_* = ', num2str(round(df.velProf.modelMeterNoAlpha.params.ustar,4))), ...
                strcat('R^2 = ', num2str(round(df.velProf.modelMeterNoAlpha.params.model.Rsquared.Ordinary,2)))};
                text(0.05, 0.7, sprintf('%s\n', params1a{:}), 'Parent', gca, 'units', 'normalized');
            end
            % plot water surface
            plot(xlim, repmat(df.FlowDepthAtCollection, 1, 2), 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.5);
            % clean up and format plot
            ylim([0 df.FlowDepthAtCollection*1.1])
            xlim([0 10])
            xlabel('velocity (m/s)')
            ylabel('distance above bed (m)')
            title(df.StationID, 'interpreter', 'none')
            box on
            set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
            
        %% grain size
        subplot(2, 2, 2); hold on
            % calculations for coloring
            [uniqueDepth, ~, uniquePlace] = unique(df.waterSamplesTable.sampleDepth);
            cmapRaw = flipud(parula( size(uniqueDepth, 1) ));
            cmap = cmapRaw(uniquePlace, :);
            % begin plotting
            semilogx(gsDistClass2num(df.gsDistBed), cumsum(df.gsDistBed{:,1}, 'omitnan'), 'k-o', 'LineWidth', 1.5, 'MarkerSize', 2)
            hold on
            for p = 1:size(df.waterSamplesTable, 1)
                semilogx(df.waterSamplesTable.gsClass{p}(2:end), cumsum(table2array(df.waterSamplesTable.gsDistNWnorm{p}(2:end,:)), 'omitnan'), ...
                    '-o', 'Color', cmap(p, :), 'LineWidth', 1.5, 'MarkerSize', 2)
            end
            plot(repmat(df.Samples(1).gsWLdata.gsWLcutoff, 2, 1), [0 100], 'k--')
            title(strrep(strtrim(df.StationID), '_', ' '));
            xlim([10 300])
            ylim([0 100])
            xlabel('grain size (\mu m)')
            ylabel('% finer')
            box on;
            set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')

        %% concentration profile
        subplot(2, 2, 3); hold on
            plot(df.waterSamplesTable.concNW, df.waterSamplesTable.sampleZ, 'ko', 'MarkerFaceColor', col3, 'MarkerSize', 5)
            plot(df.waterSamplesTable.conc, df.waterSamplesTable.sampleZ, 'ko', 'MarkerSize', 5)
            plot(df.concProf.model.Cs, df.concProf.model.Zs, 'k--', 'MarkerSize', 3, 'LineWidth', 1.5)
            plot(df.concProf.pred.CsSum, df.concProf.pred.ZsSum, 'k-', 'MarkerSize', 3, 'LineWidth', 1.5)
            text(0, 0, [df.StationID, ' ,' num2str(p)])
            text(0, 0, [num2str(df.Velocity.ustarBest)])
            ylim([0 df.FlowDepthAtCollection*1.1])
            xlim([0 25])
            xlabel('concentration (g/L)')
            ylabel('distance above bed (m)')
            title(df.StationID, 'interpreter', 'none')
            box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
        

        set(gcf, 'Pos', [50 100 900 900], 'PaperPositionMode', 'auto')
        if printOut
            tempid = df.StationID;
            tempid(df.StationID == '/') = [];
            print('-dpng', '-r300', ['./figsExport/staions2018_', tempid, '.png']);
            print('-depsc', '-painters', '-r300', ['./figsExport/staions2018_', tempid, '.eps']);
            pause(1);
        end
            
end