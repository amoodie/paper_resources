function analysis_parameterspace()


    clear all

    %% load data
    % add things to path
    source_path = genpath(fullfile('..', 'source'));
    input_path = genpath(fullfile('..', 'input_data'));
    output_data_path = genpath(fullfile('..', 'output_data'));
    plotting_path = genpath(fullfile('..', 'post_processing', 'plotting_scripts'));
    addpath(source_path, input_path, output_data_path, plotting_path);
    [CP] = color_palette();
    
    % load engineered data
    paramspace_raw = load(fullfile('paramspace_analysis.mat'));
    paramspace = paramspace_raw.paramspace;
    preavul_thresh_list = paramspace_raw.preavul_thresh_list;
    discharge_list = paramspace_raw.discharge_list;
    paramspace_setup = paramspace_raw.paramspace_setup;
    paramspace_discharge = paramspace_raw.paramspace_discharge;
    clear paramspace_raw

    %% calculations for means, stdevs, etc
    paramspaceLblow = paramspace{1}.Lblow;
    dx = paramspace{1}.cfg.dx;
    
    % fix NaN values
    paramspace_discharge(isnan(paramspace_discharge)) = 1500;
    
    % clean up any failed runs
%     failed = cell2mat(arrayfun(@(x) isempty(x{:}.avul_idx), paramspace, 'UniformOutput', false));
    failed = arrayfun(@(x) isempty(x{:}), paramspace, 'Unif', 1);
    manual_failed_idx = sub2ind(size(failed), [], []);
    failed(manual_failed_idx) = [];
    if any(failed)
        for i = find(failed)'
            paramspace{i}.avul_idx = NaN;
            paramspace{i}.avul_len = NaN;
            paramspace{i}.lobe_len = NaN;
            paramspace{i}.avul_time = NaN;
        end
    end
    
    paramspace_TA_mean = cell2mat(arrayfun((@(x)  nanmean(x{:}.avul_time(1:end) / 365)), paramspace, 'UniformOutput', false));
    paramspace_TA_std = cell2mat(arrayfun((@(x) std(x{:}.avul_time(1:end) ./ 365, 'omitnan')), paramspace, 'UniformOutput', false));
    paramspace_LA_mean = cell2mat(arrayfun((@(x)  nanmean(x{:}.avul_len(1:end) / 1000)), paramspace, 'UniformOutput', false));
    paramspace_LA_std = cell2mat(arrayfun((@(x) std(x{:}.avul_len(1:end) ./ 1000, 'omitnan')), paramspace, 'UniformOutput', false));
    paramspace_LL_mean = cell2mat(arrayfun((@(x)  nanmean(x{:}.lobe_len(1:end) / 1000)), paramspace, 'UniformOutput', false));
    paramspace_LL_std = cell2mat(arrayfun((@(x) std(x{:}.lobe_len(1:end) ./ 1000, 'omitnan')), paramspace, 'UniformOutput', false));
    
    
    %% set up the plotting vecotrs and grids
    % grid needs to be made to offset from the actual calculaiton points!
%     targets = [repmat(mapExploreSetupList', length(mapExploreTrigList), 1), ...
%         repelem(mapExploreTrigList', length(mapExploreSetupList), 1)];
    
    %% convert one of the setup threshold vectors into a table
%     demo_idx = find(preavul_thresh_list == 0.6);
%     paramspace_table = array2table([avulsion_location_list', paramspace_TA(demo_idx, :)', paramspace_TA_norm(demo_idx, :)', ...
%                              paramspace_LA(demo_idx, :)', paramspace_LL(demo_idx, :)'], ...
%                              'VariableNames', {'avulloc', 'TA', 'TA_norm', 'LA', 'LL'});

    
    
    load(fullfile('..', 'input_data', 'YRD_datatable.mat'), 'YRD_table');
    YRD_table.avul_time = [NaN; (YRD_table.date(2:end) - YRD_table.date(1:end-1)) / 365];
    % TA = 7 ± 2 yr and L A,YR = 52.5 ± 12.3 km
    
    
    %% big plot of surface data, full eps
    map_only = false;
    fig = figure();
    subplot(2, 3, 4)
        hold on
        plot_map(paramspace_TA_mean, paramspace_setup, paramspace_discharge, ...
            [7, 2, min(YRD_table.avul_time), max(YRD_table.avul_time)], map_only)
        xlabel('avul. setup thresh. (Hbf)')
        ylabel('flood discharge (m3/s)')
        h = colorbar; set(get(h, 'label'), 'string', 'avul. time (yr)'); 
        xlim([0.2 0.6])
        ylim([1500 3500])
        caxis([1 50])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
     subplot(2, 3, 5)
        hold on
        plot_map(paramspace_LA_mean, paramspace_setup, paramspace_discharge, ...
            [52.5, 12.3, min(YRD_table.avul_len), max(YRD_table.avul_len)], map_only)
        xlabel('avul. setup thresh. (Hbf)')
        ylabel('flood discharge (m3/s)')
        h = colorbar; set(get(h, 'label'), 'string', 'avul. length (km)'); 
        xlim([0.2 0.6])
        ylim([1500 3500])
        caxis([10 90])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    subplot(2, 3, 6)
        hold on
        plot_map(paramspace_LL_mean, paramspace_setup, paramspace_discharge, ...
            [0 0 0 0], map_only)
        xlabel('avul. setup thresh. (Hbf)')
        ylabel('flood discharge (m3/s)')
        h = colorbar; set(get(h, 'label'), 'string', 'lobe length (km)');
        xlim([0.2 0.6])
        ylim([1500 3500])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    set(fig, 'Pos', [200 500 1300 700]);
    set(fig, 'PaperPositionMode', 'auto')
    pause(0.1)
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs','paramspace_maps.eps'))

    
    %% big plot of surface data, empty, maps only
    map_only = true;
    fig = figure();
    subplot(2, 3, 4)
        hold on
        plot_map(paramspace_TA_mean, paramspace_setup, paramspace_discharge, ...
            [7, 2, min(YRD_table.avul_time), max(YRD_table.avul_time)], map_only)
        xlabel('avul. setup thresh. (Hbf)')
        ylabel('flood discharge (m3/s)')
        h = colorbar; set(get(h, 'label'), 'string', 'avul. time (yr)'); 
        xlim([0.2 0.6])
        ylim([1500 3500])
        caxis([1 50])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
        set(gca,'Visible','off')
     subplot(2, 3, 5)
        hold on
        plot_map(paramspace_LA_mean, paramspace_setup, paramspace_discharge, ...
            [52.5, 12.3, min(YRD_table.avul_len), max(YRD_table.avul_len)], map_only)
        xlabel('avul. setup thresh. (Hbf)')
        ylabel('flood discharge (m3/s)')
        h = colorbar; set(get(h, 'label'), 'string', 'avul. length (km)'); 
        xlim([0.2 0.6])
        ylim([1500 3500])
        caxis([10 90])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
        set(gca,'Visible','off')
    subplot(2, 3, 6)
        hold on
        plot_map(paramspace_LL_mean, paramspace_setup, paramspace_discharge, ...
            [0 0 0 0], map_only)
        xlabel('avul. setup thresh. (Hbf)')
        ylabel('flood discharge (m3/s)')
        h = colorbar; set(get(h, 'label'), 'string', 'lobe length (km)');
        xlim([0.2 0.6])
        ylim([1500 3500])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
        set(gca,'Visible','off')
    set(fig, 'Pos', [200 500 1300 700]);
    set(fig, 'PaperPositionMode', 'auto')
    drawnow
    pause(1)
    print('-dpng', '-r300', '-opengl', '../figs/paramspace_maps_only.png');
%     print('-depsc', '-painters', '-r300', fullfile('..', 'figs','paramspace_maps.eps'))

    
    %% make difference maps
    TAdiffmap = (paramspace_TA_mean - 7) ./ (7) .* 100;
    LAdiffmap = (paramspace_LA_mean - 52.5) ./ (52.5) .* 100;
    [~, TAminI] = min(min(abs(TAdiffmap(2:end,:))));
    [~, LAminI] = min(min(abs(LAdiffmap(2:end,:))));


    totdiffmap = ((paramspace_TA_mean - 7)./(7)) + ((paramspace_LA_mean - 52.500)./(52.500));
    [~, totminI] = min(min(abs(totdiffmap(2:end,:))));

    
    ta_part = ((paramspace_TA_mean - 7)./(7)).^2;
    la_part = ((paramspace_LA_mean - 52.500)./(52.500)).^2;
    ta_norm = (ta_part-min(min(ta_part)))/(max(max(ta_part))-min(min(ta_part)));
    la_norm = (la_part-min(min(la_part)))/(max(max(la_part))-min(min(la_part)));
    normdiffmap = ta_norm + la_norm ;
    
    map_only = false;
    fig = figure();
    subplot(2, 3, 4)
        hold on
        plot_map(abs(TAdiffmap), paramspace_setup, paramspace_discharge, ...
            [0, 0, 0, 0], map_only)
        xlabel('avul. setup thresh. (Hbf)')
        ylabel('flood discharge (m3/s)')
        h = colorbar; set(get(h, 'label'), 'string', 'avul. time (yr)'); 
        xlim([0.2 0.6])
        ylim([1500 3500])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
     subplot(2, 3, 5)
        hold on
        plot_map(abs(LAdiffmap), paramspace_setup, paramspace_discharge, ...
            [0 0 0 0], map_only)
        xlabel('avul. setup thresh. (Hbf)')
        ylabel('flood discharge (m3/s)')
        h = colorbar; set(get(h, 'label'), 'string', 'avul. length (km)'); 
        xlim([0.2 0.6])
        ylim([1500 3500])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    subplot(2, 3, 6)
        hold on
        plot_map(abs(normdiffmap), paramspace_setup, paramspace_discharge, ...
            [0 0 0 0], map_only)
        xlabel('avul. setup thresh. (Hbf)')
        ylabel('flood discharge (m3/s)')
        h = colorbar; set(get(h, 'label'), 'string', 'total diff map (km)');
        xlim([0.2 0.6])
        ylim([1500 3500])
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    set(fig, 'Pos', [200 500 1300 700]);
    set(fig, 'PaperPositionMode', 'auto')
    pause(0.1)
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs','paramspace_maps_diff.eps'))

        
end


function plot_map(map, X, Y, cont, map_only)

    gc = [0.85 0.85 0.85];

    % plot the main upper part
    surf(X(2:end, :), Y(2:end, :), map(2:end, :), 'EdgeColor', 'none', 'FaceColor', 'interp')
    
    % plot the lower NaN part
    surf(repmat(X(1, :), 2, 1), Y(1:2, :), repmat(map(1, :), 2, 1), 'EdgeColor', 'none', 'FaceColor', 'interp')

    if ~map_only
        % put ticks on the upper part
        plot3(X(2:end, :), Y(2:end, :), ones(size(X(2:end, :)))*max(max(map)), 'k+')

        % put ticks on the lower part
        plot3(X(1, :), mean([Y(1, :); Y(2, :)], 1), ones(size(X(2:end, :)))*max(max(map)), 'k+')

        % plot the thick line
        plot3(X(2, :), Y(2, :), ones(size(X(2:end, :)))*max(max(map)), 'k-', 'LineWidth', 2)

        % plot a contour line at the mean
        contour3(X(2:end, :), Y(2:end, :), map(2:end, :), [cont(1) cont(1)], 'linecolor', gc, 'LineWidth', 1.5, 'LineStyle', '-')
        contour3(repmat(X(1, :), 2, 1), Y(1:2, :), repmat(map(1, :), 2, 1), [cont(1) cont(1)], 'linecolor', gc, 'LineWidth', 1.5, 'LineStyle', '-')

        % plot a contour line at the stds
        contour3(X(2:end, :), Y(2:end, :), map(2:end, :), [cont(1)-cont(2) cont(1)-cont(2)], 'linecolor', gc, 'LineWidth', 1.5, 'LineStyle', '--')
        contour3(repmat(X(1, :), 2, 1), Y(1:2, :), repmat(map(1, :), 2, 1), [cont(1)-cont(2) cont(1)-cont(2)], 'linecolor', gc, 'LineWidth', 1.5, 'LineStyle', '--')
        contour3(X(2:end, :), Y(2:end, :), map(2:end, :), [cont(1)+cont(2) cont(1)+cont(2)], 'linecolor', gc, 'LineWidth', 1.5, 'LineStyle', '--')
        contour3(repmat(X(1, :), 2, 1), Y(1:2, :), repmat(map(1, :), 2, 1), [cont(1)+cont(2) cont(1)+cont(2)], 'linecolor', gc, 'LineWidth', 1.5, 'LineStyle', '--')

        % plot a contour line at the min, max
        contour3(X(2:end, :), Y(2:end, :), map(2:end, :), [cont(3) cont(3)], 'linecolor', gc, 'LineWidth', 1.5, 'LineStyle', ':')
        contour3(repmat(X(1, :), 2, 1), Y(1:2, :), repmat(map(1, :), 2, 1), [cont(3) cont(3)], 'linecolor', gc, 'LineWidth', 1.5, 'LineStyle', ':')
        contour3(X(2:end, :), Y(2:end, :), map(2:end, :), [cont(4) cont(4)], 'linecolor', gc, 'LineWidth', 1.5, 'LineStyle', ':')
        contour3(repmat(X(1, :), 2, 1), Y(1:2, :), repmat(map(1, :), 2, 1), [cont(4) cont(4)], 'linecolor', gc, 'LineWidth', 1.5, 'LineStyle', ':')
    end
end
