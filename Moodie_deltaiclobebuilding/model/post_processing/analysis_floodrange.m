function analysis_floodrange()

    clear all
    
    %% load data
    % add things to path
    source_path = genpath(fullfile('..', 'source'));
    output_data_path = genpath(fullfile('..', 'output_data'));
    plotting_path = genpath(fullfile('..', 'post_processing', 'plotting_scripts'));
    addpath(source_path, output_data_path, plotting_path);
    [CP] = color_palette();
    
%     main = load('floodrange_analysis_run1.mat');
%     front = load('floodrange_analysis_run2.mat');    
%     floodrange = main.floodrange;
%     flood_list = main.flood_list;
%     floodrange(1:2) = front.floodrange(1:2);

    main = load('floodrange_analysis05_front.mat');
    back = load('floodrange_analysis05_back1.mat');
    back2 = load('floodrange_analysis05_back2.mat');
    floodrange = main.floodrange;
    flood_list = main.flood_list;
    floodrange(end-1) = back.floodrange(1);
%     floodrange(end) = back2.floodrange(1);
    floodrange = floodrange(1:end-1);
    flood_list = flood_list(1:end-1);

    % recreate the discharge curves in order to calculate the coefficient of variation
    for i = 1:length(flood_list)
        Qw_curves(:, i) = set_Qw('engineered', flood_list(i));
        H_curves(:,i) = Hn_DW(Qw_curves(:, i), 400, 0.001, 6.4e-5, 9.81);
        CV_curves(1, i) = std(H_curves(:, i)) ./ mean(H_curves(:, i));
    end

    flood_Lblow = floodrange{1}.Lblow;
    flood_nGoodAvuls = cell2mat(cellfun((@(x) length(x.avul_time(4:end))), floodrange, 'UniformOutput', false));
    flood_TAmean = cell2mat(cellfun((@(x) nanmean(x.avul_time(4:end)/365)), floodrange, 'UniformOutput', false));
    flood_TAstd = cell2mat(cellfun((@(x) std(x.avul_time(4:end)/365, 'omitnan')), floodrange, 'UniformOutput', false));
    flood_LAmean = cell2mat(cellfun((@(x) nanmean(x.avul_len(4:end)/flood_Lblow)), floodrange, 'UniformOutput', false));
    flood_LAstd = cell2mat(cellfun((@(x) std(x.avul_len(4:end)/flood_Lblow, 'omitnan')), floodrange, 'UniformOutput', false));
    flood_LLmean = cell2mat(cellfun((@(x) nanmean(x.lobe_len(4:end)/flood_Lblow)), floodrange, 'UniformOutput', false));
    flood_LLstd = cell2mat(cellfun((@(x) std(x.lobe_len(4:end)/flood_Lblow, 'omitnan')), floodrange, 'UniformOutput', false));
    flood_table = array2table([flood_list', flood_TAmean, flood_TAstd, flood_LAmean, flood_LAstd, flood_LLmean, flood_LLstd], ...
        'VariableNames', {'flood', 'TAmean', 'TAstd', 'LAmean', 'LAstd', 'LLmean', 'LLstd'});

    for i = 1:length(floodrange)
        TA_vals(:,i) = (floodrange{i}.avul_time(4:end)/365)';
        LA_vals(:,i) = (floodrange{i}.avul_len(4:end)/flood_Lblow)';
        LL_vals(:,i) = (floodrange{i}.lobe_len(4:end)/flood_Lblow)';
    end
    jitter = ( randn(size(TA_vals))*50);
    offset = 0;
    jitterMat = repmat(flood_list, size(jitter, 1), 1) + offset + jitter;
    
    
    for i = 1:length(floodrange)
        lobelens = floodrange{i}.lobe_len(4:end)/flood_Lblow;
        zeroidx = lobelens == 0;
        lobelens(zeroidx) = NaN;
        means(i) = nanmean(lobelens);
        stds(i) = std(lobelens, 'omitnan');
    end
    flood_table.LLmeanNoZ = means';
    flood_table.LLstdNoZ = stds';

    
    fig = figure();
    subplot(2, 3, 4)
        hold on
        ymax = max(flood_table.TAmean+flood_table.TAstd) * 1.1;
        plot(jitterMat, TA_vals, 'o', 'Color', CP.g3)
        errorbar(flood_table.flood, flood_table.TAmean, flood_table.TAstd, ...
                 'Color', [0 0 0], 'LineStyle', 'none');
        plot(flood_table.flood, flood_table.TAmean, 'ok', 'MarkerFaceColor', CP.c3_s2)
%         [fig] = NaN_plotter(flood_table.flood, flood_table.TAmean, ymax, CP.c3_s2, fig);
        xlim([0 6000])
        ylim([0, ymax])
        xlabel('flood discharge (m^3/s)')
        ylabel('avulsion time (yr)')
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    subplot(2, 3, 5)
        hold on
        ymax = max(flood_table.LAmean+flood_table.LAstd) * 1.1;
        plot(jitterMat, LA_vals, 'o', 'Color', CP.g3)
        errorbar(flood_table.flood, flood_table.LAmean, flood_table.LAstd, ...
                 'Color', [0 0 0], 'LineStyle', 'none');
        plot(flood_table.flood, flood_table.LAmean, 'ok', 'MarkerFaceColor', CP.c3_s2)
%         [fig] = NaN_plotter(flood_table.flood, flood_table.LAmean, ymax, CP.c3_s2, fig);
        xlim([0 6000])
        ylim([0, ymax])
        xlabel('flood discharge (m^3/s)')
        ylabel('avulsion length (LA/Lb)')
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    subplot(2, 3, 6)
        hold on
        ymax = max(flood_table.LLmean+flood_table.LLstd) * 1.1;
        plot(jitterMat, LL_vals, 'o', 'Color', CP.g3)
        errorbar(flood_table.flood, flood_table.LLmean, flood_table.LLstd, ...
                 'Color', [0 0 0], 'LineStyle', 'none');
        plot(flood_table.flood, flood_table.LLmean, 'ok', 'MarkerFaceColor', CP.c3_s2)
%         [fig] = NaN_plotter(flood_table.flood, flood_table.LLmean, ymax, CP.c3_s2, fig);
        xlim([0 6000])
        ylim([0, ymax])
        xlabel('flood discharge (m^3/s)')
        ylabel('lobe length (LL/Lb)')
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    set(fig, 'Pos', [200 500 800 400]);
    set(fig, 'PaperPositionMode', 'auto')
    pause(0.1)
    
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs','floodrange_analysis.eps'))
        
    %% coefficient of variation plots for comparison with Chadwick 2019
    figure()
    subplot(2, 3, 2); hold on;
        ymax = max(flood_table.LAmean+flood_table.LAstd) * 1.1;
%         plot(jitterMat, LA_vals, 'o', 'Color', CP.g3)
        errorbar(CV_curves, flood_table.LAmean, flood_table.LAstd, ...
                 'Color', [0 0 0], 'LineStyle', 'none');
        plot(CV_curves, flood_table.LAmean, 'ok', 'MarkerFaceColor', CP.c3_s2)
        xlim([0 1])
        ylim([0, ymax])
        xlabel('coefficient of variation')
        ylabel('avulsion length (LA/Lb)')
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs','CV_comparison.eps'))
        
    
end

function [fig] = NaN_plotter(xvalue, yvalue, ymax, color, fig)

    nan_idx = find(isnan(yvalue));
    
    off = 100;
    for i = nan_idx
%         ylims = ylim;
        ylims = [0, ymax];
        fill([xvalue(i)-off, xvalue(i)-off, xvalue(i)+off, xvalue(i)+off], ...
             [ylims fliplr(ylims)], ...
             color, 'EdgeColor', 'none')
    end


end