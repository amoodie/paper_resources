function analysis_setuprange()

    clear all
    
    %% load data
    % add things to path
    source_path = genpath(fullfile('..', 'source'));
    output_data_path = genpath(fullfile('..', 'output_data'));
    plotting_path = genpath(fullfile('..', 'post_processing', 'plotting_scripts'));
    addpath(source_path, output_data_path, plotting_path);
    [CP] = color_palette();
    
    load('setuprange_analysis.mat')
    
    %% calculations for means, stdevs, etc
    setup_Lblow = setuprange{1}.Lblow;
    setup_TAmean = cell2mat(arrayfun((@(x) nanmean(x{:}.avul_time(4:end)/365)), setuprange, 'UniformOutput', false));
    setup_TAstd = cell2mat(arrayfun((@(x) std(x{:}.avul_time(4:end)/365, 'omitnan')), setuprange, 'UniformOutput', false));
    setup_TAmed = cell2mat(arrayfun((@(x) median(x{:}.avul_time(4:end)/365)), setuprange, 'UniformOutput', false));
    setup_LAmean = cell2mat(arrayfun((@(x) nanmean(x{:}.avul_len(4:end)/setup_Lblow)), setuprange, 'UniformOutput', false));
    setup_LAstd = cell2mat(arrayfun((@(x) std(x{:}.avul_len(4:end)/setup_Lblow, 'omitnan')), setuprange, 'UniformOutput', false));
    setup_LAmed = cell2mat(arrayfun((@(x) median(x{:}.avul_len(4:end)/setup_Lblow, 'omitnan')), setuprange, 'UniformOutput', false));
    setup_LLmean = cell2mat(arrayfun((@(x) nanmean(x{:}.lobe_len(4:end)/setup_Lblow)), setuprange, 'UniformOutput', false));
    setup_LLstd = cell2mat(arrayfun((@(x) std(x{:}.lobe_len(4:end)/setup_Lblow, 'omitnan')), setuprange, 'UniformOutput', false));
    setup_LLmed = cell2mat(arrayfun((@(x) median(x{:}.avul_len(4:end)/setup_Lblow, 'omitnan')), setuprange, 'UniformOutput', false));
    setup_table = array2table([preavul_thresh_list', setup_TAmean, setup_TAstd, setup_TAmed, setup_LAmean, setup_LAstd, setup_LAmed, setup_LLmean, setup_LLstd, setup_LLmed], ...
        'VariableNames', {'setup', 'TAmean', 'TAstd', 'TAmed', 'LAmean', 'LAstd', 'LAmed', 'LLmean', 'LLstd', 'LLmed'});
    setup_table.pidx = ~(setup_table.setup > 1.5); % plotting index for fast changes

    fig = figure();
    subplot(2, 3, 4)
        hold on
        errorbar(setup_table.setup, setup_table.TAmean, setup_table.TAstd, 'Color', [0 0 0], 'LineStyle', 'none');
        plot(setup_table.setup, setup_table.TAmean, 'ok', 'MarkerFaceColor', CP.c4_s3)
        xlim([0 1.1])
        xlabel('setup thresh. (Hbf)')
        ylabel('avulsion time (yr)')
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    subplot(2, 3, 5)
        hold on
        errorbar(setup_table.setup, setup_table.LAmean, setup_table.LAstd, 'Color', [0 0 0], 'LineStyle', 'none');
        plot(setup_table.setup, setup_table.LAmean, 'ok', 'MarkerFaceColor', CP.c4_s3)
        xlim([0 1.1])
        xlabel('setup thresh. (Hbf)')
        ylabel('avulsion length (LA/Lb)')
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    subplot(2, 3, 6)
        hold on
        errorbar(setup_table.setup, setup_table.LLmean, setup_table.LLstd, 'Color', [0 0 0], 'LineStyle', 'none');
        plot(setup_table.setup, setup_table.LLmean, 'ok', 'MarkerFaceColor', CP.c4_s3)
        xlim([0 1.1])
        ylim([0 1.5])
        xlabel('setup thresh. (Hbf)')
        ylabel('lobe length (LL/Lb)')
        box on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
    set(fig, 'Pos', [200 500 800 400]);
    set(fig, 'PaperPositionMode', 'auto')
    pause(0.1)
    
   print('-depsc', '-painters', '-r300', fullfile('..', 'figs','setuprange_analysis.eps'))
        
       

end





