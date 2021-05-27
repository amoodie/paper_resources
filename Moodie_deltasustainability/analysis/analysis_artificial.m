function analysis_artifical()

    clear all
    
    %% load data
    % add things to path
    source_path = genpath(fullfile('..', 'source'));
    input_path = genpath(fullfile('..', 'input_data'));
    output_data_path = genpath(fullfile('..', 'model', 'output_data'));
    home_path = genpath(fullfile('~', 'Documents', 'MATLAB'));
    plotting_path = genpath(fullfile('..', 'post_processing', 'plotting_scripts'));
    addpath(source_path, input_path, output_data_path, home_path, plotting_path);
    global CP
    [CP] = color_palette();
    label_interpreter = 'latex';
    printOut = true;
    print_root = './';
    bgcolor = [1 1 1]; % [0.94 0.94 0.94];
    
    % load engineered avulsion data results
    load(fullfile('engineered_single_avulsion_location_list.mat'), 'avulsion_location_list', 'n_expts')
    if false
        file_list = dir('../output_data/engineered');
        is_match = contains([{file_list.name}]', 'engineered_data_single_LA');
        match_list = [{file_list(is_match).name}]';
        avulsion_location_range = cell(length(avulsion_location_list), n_expts);
        for i = 1:size(avulsion_location_range,1)
            for j = 1:n_expts
                s = load(fullfile('..', 'output_data', 'engineered', ['engineered_data_single_LA',num2str(i),'_expt',num2str(j),'.mat']), 'data');
                s = s.data;
                d = s.d;
                [d] = extract_avulsion_information(s, d);
                [d] = extract_additional_information(s, d);
                [analysis] = convert_to_analysis(s, d);
                
                avulsion_location_range(i,j) = {analysis};

                clear data s d analysis
            end
        end
        save(fullfile('..', 'output_data', 'engineered', 'engineered_single_data_analysis.mat'), ...
            'avulsion_location_range', 'avulsion_location_list', '-v7.3')
        error('NotImplemented to flow from here')
    else
        % load it from file
        eng_raw = load('engineered/engineered_single_data_analysis.mat');
        eng.analysis = eng_raw.avulsion_location_range;
        eng.avulsion_location_list = eng_raw.avulsion_location_list;
        clear eng_raw
    end
        
    % load spinup results
    file_list = dir(fullfile('..', 'model', 'output_data', 'spinups'));
    is_match = contains([{file_list.name}]', 'engineered_spinup_analysis_single_');
    match_list = [{file_list(is_match).name}]';
    spin_raw = load(fullfile('engineered', 'spinups', match_list{1}));
    spin = spin_raw.analysis;
    spin_fn = fieldnames(spin);
    clear spin_raw
    for m = 2:length(match_list)
        spin_raw_m = load(fullfile('engineered', 'spinups', match_list{m}));
        spin_m = spin_raw_m.analysis;
        for fn = 1:length(spin_fn)
            spin.(spin_fn{fn}) = cat(2, spin.(spin_fn{fn}), spin_m.(spin_fn{fn}));
        end
        clear spin_raw_m
    end
    spin.avul_time_mean = mean(spin.avul_time/365, 'omitnan');
    spin.avul_len_mean = mean(spin.avul_len/1000, 'omitnan');
    spin.lobe_len_mean = mean(spin.lobe_len/1000, 'omitnan');
    spin.avul_time_std = std(spin.avul_time/365, [], 'omitnan');
    spin.avul_len_std = std(spin.avul_len/1000, [], 'omitnan');
    
    % load ganti experimental data
    exp_raw = readtable('delta_experiment_avulsion_lengths.csv');
    expt.tbl = exp_raw(logical(exp_raw.stable), :);
    clear exp_raw
    
    % load actual yellow river delta data
    YRD_mat = load(fullfile('YRD_datatable.mat'));
    YRD.tbl = YRD_mat.YRD_table;
    clear YRD_mat
    
    
    %% calculations for means, stdevs, etc on engineered avulsion dataset
    eng.Lblow = eng.analysis{1}.Lblow;
    dx = eng.analysis{1}.dx;
    
    % clean up any failed runs
    failed = cell2mat(arrayfun(@(x) isempty(x{:}.avul_idx), eng.analysis, 'UniformOutput', false));
    manual_failed_idx = sub2ind(size(failed), [11, 11, 11, 11, 11, 14, 14, 14, 14], [1, 2, 5, 7, 9, 2, 6, 8, 9]); % sub2ind(size(failed), [6, 6, 6, 1, 10], [6, 7, 11, 15, 11]);
    failed(manual_failed_idx) = 1;
    for i = find(failed)'
        eng.analysis{i}.avul_idx = NaN;
        eng.analysis{i}.avul_len = NaN;
        eng.analysis{i}.lobe_len = NaN;
    end
    
    % compute the means etc
    eng.TA = cell2mat(arrayfun((@(x) x{:}.avul_idx ./ 365), eng.analysis, 'UniformOutput', false));
    eng.TA_norm = eng.TA ./ spin.avul_time_mean;
    eng.LA = cell2mat(arrayfun((@(x) x{:}.avul_len ./ 1000), eng.analysis, 'UniformOutput', false));
    eng.LA_norm = eng.LA ./ spin.avul_len_mean;
    eng.LL = cell2mat(arrayfun((@(x) x{:}.lobe_len ./ 1000), eng.analysis, 'UniformOutput', false));
    eng.LL_norm = eng.LL ./ spin.lobe_len_mean;
    eng.TA_mean = nanmean(eng.TA, 2);
    eng.TA_std = std(eng.TA, [], 2, 'omitnan');
    eng.LA_mean = nanmean(eng.LA, 2);
    eng.LA_std = std(eng.LA, [], 2, 'omitnan');
    eng.LL_mean = nanmean(eng.LL, 2);
    eng.LL_std = std(eng.LL, [], 2, 'omitnan');
    
    %% convert one of the setup threshold vectors into a table
    eng.tbl = array2table([avulsion_location_list', eng.TA_mean, eng.TA_std, ...
                           eng.LA_mean, eng.LA_std, eng.LL_mean, eng.LL_std], ...
                                'VariableNames', {'avulloc', 'TA_mean', 'TA_std', 'LA_mean', 'LA_std', 'LL_mean', 'LL_std'});
    save(fullfile('..', 'model', 'output_data', 'engineered_single_analysis_summary.mat'), 'eng')
    
    %% calculations on real data
    YRD.Lb_low = 21;
    YRD.Lb_high = 54;
    YRD.Lb = 37.5;
    YRD.tbl.avul_time = [NaN; YRD.tbl.date(2:end) - YRD.tbl.date(1:end-1)] ./ 365;
    YRD.tbl.next_avul_time = [YRD.tbl.date(2:end) - YRD.tbl.date(1:end-1); NaN] ./ 365;
    YRD.tbl.next_avul_time_norm = YRD.tbl.next_avul_time ./ nanmean(YRD.tbl.avul_time);
    YRD.tbl.next_avul_time_norm_low = YRD.tbl.next_avul_time ./ (nanmean(YRD.tbl.avul_time)+std(YRD.tbl.avul_time));
    YRD.tbl.next_avul_time_norm_high = YRD.tbl.next_avul_time ./ (nanmean(YRD.tbl.avul_time)-std(YRD.tbl.avul_time));
    
    YRD.tbl.avul_len_norm = YRD.tbl.avul_len ./ (nanmean(YRD.tbl.avul_len));
    YRD.tbl.avul_len_norm_Lb = YRD.tbl.avul_len ./ YRD.Lb;
    YRD.tbl.avul_len_norm_low = YRD.tbl.avul_len ./ (nanmean(YRD.tbl.avul_len) + std(YRD.tbl.avul_len));
    YRD.tbl.avul_len_norm_high = YRD.tbl.avul_len ./ (nanmean(YRD.tbl.avul_len) - std(YRD.tbl.avul_len));
    
    YRD.tbl.next_avul_len = [YRD.tbl.avul_len(2:end); NaN];
    YRD.tbl.next_avul_len_norm = YRD.tbl.next_avul_len ./ (nanmean(YRD.tbl.avul_len));
    YRD.tbl.avul_len_norm_low = YRD.tbl.next_avul_len ./ (nanmean(YRD.tbl.avul_len) + std(YRD.tbl.avul_len));
    YRD.tbl.avul_len_norm_high = YRD.tbl.next_avul_len ./ (nanmean(YRD.tbl.avul_len) - std(YRD.tbl.avul_len));
    
    
    %% calculations on experimental delta
    expt.Lb = 2900;
    expt.tbl.next_avul_time = [expt.tbl.avul_time(2:end); NaN];
    expt.tbl.avul_len_norm_Lb = expt.tbl.avul_len ./ expt.Lb;
    expt.tbl.next_avul_len = [expt.tbl.avul_len(2:end); NaN];
    expt.tbl.next_avul_len_norm = expt.tbl.next_avul_len ./ mean(expt.tbl.avul_len);
    % nondimensionalizing time in the experiment, use Chadwick scaling param
    %          40min   36g   1000mm3   1     60min     1        1           15min   60g   1000mm3   1      60min     1        1 
    % ----  = (----- * --- * ------- * --- * -----) * (------ * --- )   +  (----- * --- * ------- * ---  * -----) * (------ * --- )
    %  hr      60min   min   1.30g     70mm  1hr       2900mm   9.5mm       60min   min   1.30g     70mm   1hr       2900mm   13mm
    %          
    %        frac run  Qs    density   width  time     Lb      Hc
    expt.per_hour = ((40/60*36*1000/1.3/70*60) * (1/2900/9.5)) + ((15/60*60*1000/1.3/70*60) * (1/2900/13)); % scaling factor for run time (in per hour run time
    expt.tbl.next_avul_time_normTc = expt.tbl.next_avul_time / expt.per_hour;
    expt.tbl.next_avul_time_norm = expt.tbl.next_avul_time ./ mean(expt.tbl.avul_time);
    
    %% scatter plots
    mksize = 40;
    pcmap = summer(size(eng.tbl,1));
    pcmap = repmat(pcmap(12,:), size(eng.tbl,1), 1);
    if false
        fig = figure();
        subplot(2, 3, 4)
            hold on
            errorbar(eng.tbl.avulloc, eng.tbl.TA_mean ./ spin.avul_time_mean, eng.tbl.TA_std ./ spin.avul_time_mean, 'Color', [0 0 0], 'LineStyle', 'none')
            scatter(eng.tbl.avulloc, eng.tbl.TA_mean ./ spin.avul_time_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
            xlabel('artificial avul. length ($L_A/L_b$)', 'interpreter', label_interpreter)
            ylabel('next avul. time ($T_A/T_{A0}$)', 'interpreter', label_interpreter)
            xlim([0 1.7])
            ylim([0 max(eng.tbl.TA_mean ./ spin.avul_time_mean)*1.1])
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
        subplot(2, 3, 5)
            hold on
            errorbar(eng.tbl.avulloc, eng.tbl.LA_mean ./ spin.avul_len_mean, eng.tbl.LA_std ./ spin.avul_len_mean, 'Color', [0 0 0], 'LineStyle', 'none')
            scatter(eng.tbl.avulloc, eng.tbl.LA_mean ./ spin.avul_len_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
            xlim([0 1.7])
            ylim([0 max(eng.tbl.LA_mean ./ spin.avul_len_mean)*1.1])
            xlabel('artificial avul. length ($L_A/L_b$)', 'interpreter', label_interpreter)
            ylabel('next avul. length ($L_A/L_{A0}$)', 'interpreter', label_interpreter)
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
        subplot(2, 3, 6)
            hold on
            errorbar(eng.tbl.avulloc, eng.tbl.LL_mean ./ spin.Lblow/1000, eng.tbl.LL_std ./ (spin.Lblow/1000), 'Color', [0 0 0], 'LineStyle', 'none')
            scatter(eng.tbl.avulloc, eng.tbl.LL_mean ./ (spin.Lblow/1000), mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
            xlim([0 1.7])
            ylim([0 max(eng.tbl.LL_mean ./ spin.lobe_len_mean)*1.1])
            xlabel('artificial avul. length ($L_A/L_b$)', 'interpreter', label_interpreter)
            ylabel('next lobe length ($L_L/L_{L0}$)', 'interpreter', label_interpreter)
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
        set(fig, 'Pos', [200 500 800 400]);
        set(fig, 'PaperPositionMode', 'auto')
        pause(0.1)
        print('-depsc', '-painters', '-r300', fullfile('..', 'figs','engineered_analysis_nondims.eps'))

        fig = figure();
        subplot(2, 3, 4)
            hold on
            plot([0 1.7], [spin.avul_time_mean, spin.avul_time_mean], '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.5)
            errorbar(eng.tbl.avulloc, eng.tbl.TA_mean, eng.tbl.TA_std, 'Color', [0 0 0], 'LineStyle', 'none')
            scatter(eng.tbl.avulloc, eng.tbl.TA_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
            xlabel('artificial avul. length ($L_A/L_b$)', 'interpreter', label_interpreter)
            ylabel('next avul. time ($T_A$)', 'interpreter', label_interpreter)
            xlim([0 1.7])
            ylim([0 max(eng.tbl.TA_mean)*1.1])
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
        subplot(2, 3, 5)
            hold on
            plot([0 1.7], [spin.avul_len_mean, spin.avul_len_mean], '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.5)
            errorbar(eng.tbl.avulloc, eng.tbl.LA_mean, eng.tbl.LA_std, 'Color', [0 0 0], 'LineStyle', 'none')
            scatter(eng.tbl.avulloc, eng.tbl.LA_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
            xlim([0 1.7])
            ylim([0 max(eng.tbl.LA_mean)*1.1])
            xlabel('artificial avul. length ($L_A/L_b$)', 'interpreter', label_interpreter)
            ylabel('next avul. length ($L_A$)', 'interpreter', label_interpreter)
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
        subplot(2, 3, 6)
            hold on
            plot([0 1.7], [spin.lobe_len_mean, spin.lobe_len_mean], '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.5)
            errorbar(eng.tbl.avulloc, eng.tbl.LL_mean, eng.tbl.LL_std, 'Color', [0 0 0], 'LineStyle', 'none')
            scatter(eng.tbl.avulloc, eng.tbl.LL_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
            xlim([0 1.7])
            ylim([0 max(eng.tbl.LL_mean)*1.1])
            xlabel('artificial avul. length ($L_A/L_b$)', 'interpreter', label_interpreter)
            ylabel('next lobe length ($L_L$)', 'interpreter', label_interpreter)
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
        set(fig, 'Pos', [200 500 800 400]);
        set(fig, 'PaperPositionMode', 'auto')
        pause(0.1)
        print('-depsc', '-painters', '-r300', fullfile('..', 'figs','engineered_analysis_dims.eps'))

        %% scatter plot with full data sources
        fig = figure();
        subplot(2, 3, 4)
            hold on
            plot(eng.tbl.avulloc, eng.tbl.TA_mean ./ spin.avul_time_mean, 'ok', 'MarkerFaceColor', CP.c1_s1)
            plot(expt.tbl.avul_len_norm_Lb, expt.tbl.next_avul_time_norm, 'ok', 'MarkerFaceColor', CP.c3_s1)
            plot(YRD.tbl.avul_len_norm_Lb, YRD.tbl.next_avul_time_norm, 'ok', 'MarkerFaceColor', CP.c2_s1)
            xlabel('artificial avul. length ($L_A/L_b$)', 'interpreter', label_interpreter)
            ylabel('next avul. time ($T_A/T_{A0}$)', 'interpreter', label_interpreter)
            xlim([0 1.7])
            legend('model', 'exper.', 'YRD')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
        subplot(2, 3, 5)
            hold on
            plot(eng.tbl.avulloc, eng.tbl.LA_mean ./ spin.avul_len_mean, 'ok', 'MarkerFaceColor', CP.c1_s1)
            plot(expt.tbl.avul_len_norm_Lb, expt.tbl.next_avul_len_norm, 'ok', 'MarkerFaceColor', CP.c3_s1)
            plot(YRD.tbl.avul_len_norm_Lb, YRD.tbl.next_avul_len_norm, 'ok', 'MarkerFaceColor', CP.c2_s1)
            xlim([0 1.7])
            xlabel('artificial avul. length ($L_A/L_b$)', 'interpreter', label_interpreter)
            ylabel('next avul. length ($L_A/L_{A0}$)', 'interpreter', label_interpreter)
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
        subplot(2, 3, 6)
            hold on
            plot(eng.tbl.avulloc, eng.tbl.LL_mean ./ spin.lobe_len_mean, 'ok', 'MarkerFaceColor', CP.c1_s1)
            xlim([0 1.7])
            xlabel('artificial avul. length ($L_A/L_b$)', 'interpreter', label_interpreter)
            ylabel('next lobe length ($L_L/L_{L0}$)', 'interpreter', label_interpreter)
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
        set(fig, 'Pos', [200 500 800 400]);
        set(fig, 'PaperPositionMode', 'auto')
        pause(0.1)
        print('-depsc', '-painters', '-r300', fullfile('..', 'figs','engineered_analysis_alldata.eps')) 
        
        % dimensional selected plot
        fig = figure();
        subplot(1, 2, 1)
            hold on
            fill([0 0 1.7 1.7], [spin.avul_time_mean-spin.avul_time_std, spin.avul_time_mean+spin.avul_time_std, ...
                                 spin.avul_time_mean+spin.avul_time_std, spin.avul_time_mean-spin.avul_time_std], [0.9, 0.9, 0.9], 'EdgeColor', 'none')
            plot([0 1.7], [spin.avul_time_mean, spin.avul_time_mean], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
            text(0.1, spin.avul_time_mean+1.25, '$T_{A,0}$', 'Color', [0.5, 0.5, 0.5], 'FontSize', 10)
            errorbar(eng.tbl.avulloc, eng.tbl.TA_mean, eng.tbl.TA_std, 'Color', [0 0 0], 'LineStyle', 'none')
            scatter(eng.tbl.avulloc, eng.tbl.TA_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
            text(0.85, 0.1, 'a', 'units', 'normalized', 'FontSize', 12)
            ylabel('next avul. time ${T_A}''$ (yr)', 'interpreter', label_interpreter)
            xlim([0 1.7])
            ylim([0 max(eng.tbl.TA_mean)*1.1])
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
            axis square
        subplot(1, 2, 2)
            hold on
            fill([0 0 1.7 1.7], [spin.avul_len_mean-spin.avul_len_std, spin.avul_len_mean+spin.avul_len_std, ...
                                 spin.avul_len_mean+spin.avul_len_std, spin.avul_len_mean-spin.avul_len_std], [0.9, 0.9, 0.9], 'EdgeColor', 'none')
            plot([0 1.7], [spin.avul_len_mean, spin.avul_len_mean], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
            text(1.4, spin.avul_len_mean-5, '$L_{A,0}$', 'Color', [0.6, 0.6, 0.6], 'FontSize', 10)
            errorbar(eng.tbl.avulloc, eng.tbl.LA_mean, eng.tbl.LA_std, 'Color', [0 0 0], 'LineStyle', 'none')
            scatter(eng.tbl.avulloc, eng.tbl.LA_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
            text(0.85, 0.1, 'b', 'units', 'normalized', 'FontSize', 12)
            xlim([0 1.7])
            ylim([0 max(eng.tbl.LA_mean)*1.3])
            xl = xlabel('artificial avulsion length (${L_A}^* = L_A/L_b$)', 'interpreter', label_interpreter);
            ylabel('next avul. length ${L_A}''$ (km)', 'interpreter', label_interpreter)
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            set(gca, 'Layer', 'top')
            axis square
            xl.Position(1) = xl.Position(1) - 1.25;
            xl.Position(2) = xl.Position(2) - 3.5;
        set(fig, 'Pos', [200 500 570 260]); %     set(fig, 'Pos', [200 500 500 400]);
        set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        pause(0.1)
        if printOut
            print('-dpng', '-r300', [print_root, 'engineered_analysis.png']);
            set(gcf, 'InvertHardcopy', 'on')
            print('-depsc', '-painters', '-r300', [print_root, 'engineered_analysis.eps']);
            pause(1);
        end
    end
    
    %% selected scatter plots
    fig = figure();
    subplot(1, 2, 1)
        hold on
        fill([0 0 1.7 1.7], [(spin.avul_time_mean-spin.avul_time_std)/spin.avul_time_mean, (spin.avul_time_mean+spin.avul_time_std)/spin.avul_time_mean, ...
                             (spin.avul_time_mean+spin.avul_time_std)/spin.avul_time_mean, (spin.avul_time_mean-spin.avul_time_std)/spin.avul_time_mean], [0.9, 0.9, 0.9], 'EdgeColor', 'none')
        plot([0 1.7], [1, 1], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        plot((0:0.05:1.9), (0:0.05:1.9).^0.5, '-', 'Color', CP.c2_s3, 'LineWidth', 3)
        text(0.5, 0.45, '$\sqrt{{L_D}^*}$', 'rotation', 30, 'Color', CP.c2_s1, 'fontweight', 'bold', 'fontsize', 12)
        errorbar(eng.tbl.avulloc, eng.tbl.TA_mean./spin.avul_time_mean, (eng.tbl.TA_std)./spin.avul_time_mean, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)
        scatter(eng.tbl.avulloc, eng.tbl.TA_mean./spin.avul_time_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
        text(0.85, 0.1, 'a', 'units', 'normalized', 'FontSize', 12)
        ylabel('time to subsequent avulsion (${{T_A}^*}$)', 'interpreter', label_interpreter)
        yticks(0:0.2:1.6)
        xlim([0 1.7])
        ylim([0 1.6])
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        set(gca, 'Layer', 'top')
        axis square
    subplot(1, 2, 2)
        hold on
        fill([0 0 1.7 1.7], [(spin.avul_len_mean-spin.avul_len_std)/spin.avul_len_mean, (spin.avul_len_mean+spin.avul_len_std)/spin.avul_len_mean, ...
                             (spin.avul_len_mean+spin.avul_len_std)/spin.avul_len_mean, (spin.avul_len_mean-spin.avul_len_std)/spin.avul_len_mean], [0.9, 0.9, 0.9], 'EdgeColor', 'none')
        plot([0 1.7], [1 1], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        errorbar(eng.tbl.avulloc, eng.tbl.LA_mean/spin.avul_len_mean, (eng.tbl.LA_std)/spin.avul_len_mean, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)
        scatter(eng.tbl.avulloc, eng.tbl.LA_mean/spin.avul_len_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
        text(0.85, 0.1, 'b', 'units', 'normalized', 'FontSize', 12)
        xlim([0 1.7])
        ylim([0 1.6])
        xl = xlabel('diversion length (${L_D}^* = L_D/L_b$)', 'interpreter', label_interpreter);
        ylabel('subsequent avulsion length (${{L_A}^*}$)', 'interpreter', label_interpreter)
        yticks(0:0.2:1.6)
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        set(gca, 'Layer', 'top')
        axis square
        xl.Position(1) = xl.Position(1) - 1.25;
        xl.Position(2) = xl.Position(2) - 0.075;
    set(fig, 'Pos', [200 500 570 270]); %     set(fig, 'Pos', [200 500 500 400]);
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    pause(0.1)
    if printOut
        print('-dpng', '-r300', [print_root, 'engineered_analysis.png']);
        set(gcf, 'InvertHardcopy', 'on')
        print('-depsc', '-opengl', '-r600', [print_root, 'engineered_analysis.eps']);
        pause(1);
    end
    
    
    %% movie of the full spin-up
    if false
        movFig = figure();
        movie_colored(movFig, fullfile('engineered', 'spinups', ['engineered_spinup_data_single.mat']), 300, fullfile('engineered', 'spinup'));
    end
    
    
    %% computations on wave of bed erosion for each case
    if false
        series = cell(size(eng.analysis));
        for p = 1:numel(series)
            [i, j] = ind2sub(size(eng.analysis), p);
            s = load(fullfile('..', 'output_data', 'engineered', ['engineered_data_single_LA',num2str(i),'_expt',num2str(j),'.mat']), 'data');
            
            avulloc_idx = s.data.avul(end); %where the cycle-ending avulsion occurs
            taper_idx = s.data.idx_upstream_avulsion; % upstream end of taper, above where the cycle-starting diversion occurs
            eta_avulloc = s.data.eta(avulloc_idx,:);
            Hnbf = s.data.d.Hnbf;
            superelev = 1 - ((s.data.zed(avulloc_idx,:) - s.data.eta(avulloc_idx,:))./Hnbf);
            superelev_1yr = superelev(365);
            norm_superelev = superelev ./ superelev(end-1);
            [~, superelev_min_idx] = min(superelev(1:end-1));
            erosion_depth_grid = s.data.eta(1:taper_idx+11, 1) - s.data.eta(1:taper_idx+11, 2:end-1);
            [v_met, i_met] = max(erosion_depth_grid, [], 1);
            elevation_change_avulloc = ( s.data.eta(avulloc_idx,end-1)-s.data.eta(avulloc_idx,:) ) ./ Hnbf;
            [elevation_change_avulloc_max, elevation_change_avulloc_max_idx] = max(elevation_change_avulloc(1:end-1));
            % i_met =  index of max erosion at all times
            % v_met =  value of max erosion at all times
            [~, t_me] = max(v_met); % time of max erosion
            ero_depth_max = s.data.eta(i_met(t_me), 1) - s.data.eta(i_met(t_me), t_me);
            erosion_wave_ratio = ( erosion_depth_grid ) ./ ero_depth_max;
            [I,J] = (find(erosion_wave_ratio>(1/exp(1)))); % erosion_wave_contour  
            erosion_wave_upstream_index = min(I); % most upstream index
            erosion_wave_upstream_index_idxs = J(I==erosion_wave_upstream_index);
            erosion_wave_duration = max(0, max(erosion_wave_upstream_index_idxs) / 365);
            erosion_wave_len = (dx * (i_met(t_me) - erosion_wave_upstream_index) ) / s.data.d.Lblow;
            if erosion_wave_len <= 0
                erosion_wave_len = NaN;
            end
            series(i,j) = {struct('eta_avulloc', eta_avulloc, 'Hnbf', Hnbf, ...
                           'superelev', superelev, 'superelev_1yr', superelev_1yr, ...
                           'superelev_min_idx', superelev_min_idx, 'erosion_wave_len', erosion_wave_len, ...
                           'erosion_wave_depth', ero_depth_max, 'erosion_wave_dur', erosion_wave_duration, ...
                           'norm_superelev', norm_superelev, 'elevation_change_avulloc_max', elevation_change_avulloc_max)};
            clear data s d analysis
        end
        save(fullfile('..', 'output_data', 'engineered', 'engineered_single_data_erosionseries.mat'), ...
            'series', '-v7.3')
    else
        % load it from file
        ser_raw = load('engineered/engineered_single_data_erosionseries.mat');
        series = ser_raw.series;
        clear ser_raw
    end
    
    %% three long profile figure
    prows = [1, 10, 15];
    pcols = [3, 6, 4];
    pmatch = sub2ind(size(eng.analysis), prows, pcols);
    if true
        pfig = figure();
        set(pfig, 'Pos', [2000 500 800 350]);
        
        for p = 1:length(pmatch)
            [i, j] = ind2sub(size(eng.analysis), pmatch(p));
            s = load(fullfile('..', 'output_data', 'engineered', ['engineered_data_single_LA',num2str(i),'_expt',num2str(j),'.mat']), 'data');

            data_path.s = s.data;
            data_path.d = s.data.d;
            axidx = find(prows == i); % should only ever be one element (or fix prows def)
            data_path.avulloc = avulsion_location_list(prows(axidx));    
            [pfig] = plot_avulsion_cycle_records(data_path, pfig, axidx);
        end
        set(pfig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        pause(0.1)
        if true
            print('-dpng', '-r300', [print_root, 'engineered_scour_record.png']);
            set(fig, 'InvertHardcopy', 'on')
            print('-depsc', '-opengl', '-r300', [print_root, 'engineered_scour_record.eps']);
            pause(1);
        end
        % close(pfig)
    end
    
    %% fig of various erosion properties
    erosion_wave_len_mean = nanmean(cellfun(@(x) x.erosion_wave_len, series),2);
    erosion_wave_len_std = std(cellfun(@(x) x.erosion_wave_len, series),[],2);
    erosion_wave_depth_mean = nanmean(cellfun(@(x) x.erosion_wave_depth, series),2);
    erosion_wave_depth_std = std(cellfun(@(x) x.erosion_wave_depth, series),[],2);
    erosion_wave_dur_mean = nanmean(cellfun(@(x) x.erosion_wave_dur, series),2);
    erosion_wave_dur_std = std(cellfun(@(x) x.erosion_wave_dur, series),[],2);
    erosion_wave_Hnbf = mean(cellfun(@(x) x.Hnbf, series),2);
    avul_cycle_dur_mean = nanmean(cellfun(@(x) length(x.superelev)./365, series),2);
    
    erosion_wave_avulloc_max_mean = mean(cellfun(@(x) x.elevation_change_avulloc_max, series),2);
    erosion_wave_avulloc_max_std = std(cellfun(@(x) x.elevation_change_avulloc_max, series),[],2);
    
    fig = figure();
    subplot(2, 2, 1); hold on;
        plot([0 1.7], [0, 1.7], '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.2)
        text(0.2, 0.35, '1:1 line', 'interpreter', label_interpreter, 'Color', [0.4, 0.4, 0.4], 'rotation', 42)
        errorbar(eng.tbl.avulloc, erosion_wave_len_mean, erosion_wave_len_std, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)
        scatter(eng.tbl.avulloc, erosion_wave_len_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
        text(0.85, 0.1, 'a', 'units', 'normalized', 'FontSize', 12)
        ylabel({'scour length', '(${L_s}^* =L_s/L_b$)'}, 'interpreter', label_interpreter);
        yticks(0:0.2:1.4)
        xlim([0 1.7])
        ylim([0 1])
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        set(gca, 'Layer', 'top')
        axis square
    subplot(2, 2, 2); hold on;
        errorbar(eng.tbl.avulloc, erosion_wave_depth_mean./erosion_wave_Hnbf, erosion_wave_depth_std./erosion_wave_Hnbf, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)
        scatter(eng.tbl.avulloc, erosion_wave_depth_mean./erosion_wave_Hnbf, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
        text(0.85, 0.1, 'a', 'units', 'normalized', 'FontSize', 12)
        ylabel({'scour depth', '(${H_s}^* = H_s/H_{bf}$)'}, 'interpreter', label_interpreter);
        yticks(0:0.2:1.4)
        xlim([0 1.7])
        ylim([0 0.5])
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        set(gca, 'Layer', 'top')
        axis square
     subplot(2, 2, 3); hold on;
        errorbar(eng.tbl.avulloc, erosion_wave_dur_mean./avul_cycle_dur_mean, erosion_wave_dur_std./avul_cycle_dur_mean, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)
        scatter(eng.tbl.avulloc, erosion_wave_dur_mean./avul_cycle_dur_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
        text(0.85, 0.1, 'a', 'units', 'normalized', 'FontSize', 12)
        ylabel({'scour duration', '(${T_s}^* = T_s/T_{A}$)'}, 'interpreter', label_interpreter);
        yticks(0:0.2:1.4)
        xlim([0 1.7])
        ylim([0 0.5])
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        set(gca, 'Layer', 'top')
        axis square
    subplot(2, 2, 4); hold on;
        errorbar(eng.tbl.avulloc, erosion_wave_avulloc_max_mean, erosion_wave_avulloc_max_std, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)
        scatter(eng.tbl.avulloc, erosion_wave_avulloc_max_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
        text(0.85, 0.1, 'a', 'units', 'normalized', 'FontSize', 12)
        ylabel({'scour depth at ${L_A}{''}$', '(${H_s}^* = H_s/H_{bf}$)'}, 'interpreter', label_interpreter);
        yticks(0:0.2:1.4)
        xlim([0 1.7])
        ylim([0 0.5])
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        set(gca, 'Layer', 'top')
        axis square
        set(fig, 'Pos', [200 500 570 560]); %     set(fig, 'Pos', [200 500 500 400]);
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    pause(0.1)
    if printOut
        print('-dpng', '-r300', [print_root, 'engineered_scour_supp.png']);
        set(gcf, 'InvertHardcopy', 'on')
        print('-depsc', '-painters', '-r300', [print_root, 'engineered_scour_supp.eps']);
        pause(1);
    end
    
        
    %% SCOUR AND TWO LONG PROF selected plot
    fig = figure();
    % add the long profile plots for two selected indicies
    prows = [3, 12];
    pcols = [4, 6];
    pmatch = sub2ind(size(eng.analysis), prows, pcols);
    for p = 1:length(pmatch)
        [i, j] = ind2sub(size(eng.analysis), pmatch(p));
        s = load(fullfile('..', 'output_data', 'engineered', ['engineered_data_single_LA',num2str(i),'_expt',num2str(j),'.mat']), 'data');

        data_path.s = s.data;
        data_path.d = s.data.d;
        axidx = find(prows == i);
        data_path.avulloc = avulsion_location_list(prows(axidx));  
        [fig] = plot_single_avulsion_cycle_long_profile(data_path, fig, p);
    end
    subplot(2, 2, 3)
        hold on
        plot([0 1.7], [0, 1.7], '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.2)
        text(0.6, 0.7, '1:1 line', 'interpreter', label_interpreter, 'Color', [0.4, 0.4, 0.4], 'rotation', 60)
        errorbar(eng.tbl.avulloc, erosion_wave_len_mean, erosion_wave_len_std, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)
        scatter(eng.tbl.avulloc, erosion_wave_len_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
        text(0.1, 0.85, 'c', 'units', 'normalized', 'FontSize', 12)
        ylabel({'scour length (${L_s}^*$)'}, 'interpreter', label_interpreter);
        yticks(0:0.2:1.4)
        xlim([0 1.7])
        ylim([0 1])
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        set(gca, 'Layer', 'top')
        axis square
    subplot(2, 2, 4)
        hold on
        fill([repmat(1-spin.dx(1).*11./spin.Lblow(1), 2, 1); 1.7; 1.7], [0 2 2 0], ...
              [0.9, 0.9, 0.9], 'EdgeColor', 'none')
        text(((1-spin.dx(1).*11./spin.Lblow(1))+1.65)/2, 0.8, {'artificially', 'lowered'}, 'interpreter', label_interpreter, 'Color', [0, 0, 0], 'HorizontalAlignment', 'center')
        plot([0 1.7], [0, 1.7], '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.2)
        text(0.30, 0.42, '1:1 line', 'interpreter', label_interpreter, 'Color', [0.4, 0.4, 0.4], 'rotation', 56)
        errorbar(eng.tbl.avulloc, erosion_wave_avulloc_max_mean, erosion_wave_avulloc_max_std, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)
        scatter(eng.tbl.avulloc, erosion_wave_avulloc_max_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
        text(0.1, 0.85, 'd', 'units', 'normalized', 'FontSize', 12)
        xl1 = xlabel('diversion length (${L_D}^* = L_D/L_b$)', 'interpreter', label_interpreter);
        ylabel({'scour depth (${H_s}^*$)'}, 'interpreter', label_interpreter);
        yticks(0:0.2:1)
        xlim([0 1.7])
        ylim([0 1])
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        set(gca, 'Layer', 'top')
        axis square
        xl1.Position(1) = xl1.Position(1) - 1.58;
        xl1.Position(2) = xl1.Position(2) + 0.11;
    
    set(fig, 'Pos', [200 500 570 552]);
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    pause(0.1)
    if printOut
        print('-dpng', '-r300', [print_root, 'engineered_scour_waves.png']);
        set(gcf, 'InvertHardcopy', 'on')
        print('-depsc', '-painters', '-r600', [print_root, 'engineered_scour_waves.eps']);
        pause(1);
        export_fig([print_root, 'engineered_scour_waves.pdf'])
        pause(1);
    end
        
    %% figure of erosional wave
    fig = figure(); hold on;
        pmap = summer(size(eng.tbl,1));
        plot([0, 1], [0.5 0.5], '--', 'LineWidth', 1.2, 'Color', [0.6 0.6 0.6])
        text(0.53, 0.505, 'superelevation = 0.5', 'interpreter', label_interpreter, 'Color', [0.4, 0.4, 0.4])
        prows = [1, 3, 5, 7, 9];% , 11, 13];
        pcols = [2, 8, 6, 6, 2];%, 2, 1];
        plot_list = sub2ind(size(series), prows, pcols);
        for p = 1:length(plot_list)
            pp = plot_list(p);
            [lne(p)] = plot(linspace(0, 1, length(series{pp}.eta_avulloc)-1), series{pp}.norm_superelev(1:end-1), 'color', pmap(prows(p),:), 'LineWidth', 2);
        end
        ylim([0.30, 0.52])
        xlabel('normalized time ($t/{T_A}$)', 'interpreter', label_interpreter)
        ylabel('superelevation', 'interpreter', label_interpreter)
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
        axis square
        sppos = get(gca, 'position');
        subax1 = axes('Position', [sppos(1)+sppos(3)-0.34, sppos(2)+0.14, 0.25, 0.29]); hold on;
            plot([0 1.7], [0, 1.7], '--', 'Color', [0.4, 0.4, 0.4])
            errorbar(eng.tbl.avulloc, erosion_wave_len_mean, erosion_wave_len_std, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)
            scatter(eng.tbl.avulloc, erosion_wave_len_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
            text(0.2, 0.05, '1:1 line', 'interpreter', label_interpreter, 'Color', [0.4, 0.4, 0.4], 'rotation', 58)
            ylabel(subax1, {'scour len.', '(${L_s}^* =L_s/L_b$)'}, 'interpreter', label_interpreter);
            xlabel(subax1, '${L_A}^*$', 'interpreter', label_interpreter);
            xlim([0 1.7])
            ylim([0 1.0])
            axis square
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            box on
            set(gca, 'Layer', 'top')
            drawnow; pause(0.1);
    set(fig, 'Pos', [200 500 460 320]);
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    pause(0.1)
    if printOut
        print('-dpng', '-r300', [print_root, 'engineered_scour_depth.png']);
        set(fig, 'InvertHardcopy', 'on')
        print('-depsc', '-painters', '-r300', [print_root, 'engineered_scour_depth.eps']);
        pause(1);
    end
    
    
    %% do some things in loops for movies
    nframes = 75;
    for p = 1:length(pmatch)
        
        [i, j] = ind2sub(size(eng.analysis), pmatch(p));
        s = load(fullfile('..', 'output_data', 'engineered', ['engineered_data_single_LA',num2str(i),'_expt',num2str(j),'.mat']), 'data');
        
        data_path.s = s.data;
        data_path.d = s.data.d;
        
        if false
            movFig = figure();
            movie_colored(movFig, data_path, nframes, fullfile('engineered', num2str(p)));
        end
        
        dur_sim(p) = length(s.data.Qw);
        
    end
    
    t_l = 15; % time of longest, 15 seconds
    fps_mov = nframes ./ ((dur_sim / max(dur_sim)) * t_l)
    
    
    %% spinup movie
    if true
        %s = load(fullfile('..', 'output_data', 'engineered', 'spinups', 'engineered_spinup_data_single_hist_3.mat'));
        
        %data_path.s = s.data;
        %data_path.d = s.data.d;
        
        if true
            movFig = figure();
            pathstr = fullfile('..', 'output_data', 'engineered', 'spinups', 'engineered_spinup_data_single_hist_3.mat');
            movie_colored(movFig, pathstr, 50, fullfile('engineered', 'spinup'));
        end
    end
    
end


function [pfig] = plot_single_avulsion_cycle_long_profile(data_struct, pfig, axidx)

    s = data_struct.s;
    d = s.d;
    global CP

    x = d.x;
    dx = d.dx;
    Hnbf = d.Hnbf;
    L = max(d.x);
    nx = size(x, 2);
    
    % params
    upperlim = size(s.eta, 2)-1;
    nt = upperlim ./ 365;
    nprint = 10;
    
    % setup
    pidx = [1:floor(upperlim/nprint):upperlim, upperlim];
    colmap = parula(length(pidx));
    
    % colors
    Cwater1 = [195, 225, 255] ./ 255;
    Cwater2 = [35, 174, 216] ./ 255;
    Cland1 = [147, 236, 140] ./ 255;
    Cland2 = [123, 208, 117] ./ 255; 
    Csubstr1 = [196, 181, 132] ./ 255;
    Csubstr2 = [100, 63, 15] ./ 255;
    Cdgray = [0.5 0.5 0.5];
    
    % initial values
    eta_i = s.eta(:, 1);
    zed_i = s.zed(:, 1);
    cfg = yellowriver_configuration();
    [eta_base, ~, ~] = setup_initial_channel('new', floor((nx-1)/2), Hnbf, cfg.S0, nx-1, dx, ...
                                             repmat(cfg.Bc0, 1, nx), cfg.H0, cfg.Cf, ...
                                             cfg.Qwbf, cfg.eta_basin);
    
    % datums
    baselinevect = -20 * ones(nx, 1);
    gammapexxy = [cosd(d.gamm/2)*(d.gammapex), sind(d.gamm/2)*(d.gammapex)];
    
    x_norm = x - s.rad(1); % normalization x vector for rad_i for long profiles
    max_x = ((max(s.mou_idx)*dx)- s.rad(1))*1.5; % max x limit for long profiles
    
    [plan_i] = make_planform(x, d.gamm, d.gammapex, s.rad(1), false);
    xy_norm = [plan_i(end, 1), plan_i(1, 2)]; % normalization x and y scalar for planform plots
    [plan_i_norm] = plan_i - xy_norm;
    
    %% calculations
    % calculate the aggradation metrics
    agg = s.eta(:, 2:end) - s.eta(:, 1:end-1);
    rad_idx_i = get_idx(x, s.rad(1));
    avul_pt = s.avul(end);

    %% long prof
    sp1 = subplot(2, 2, axidx); % long profile
    hold on

    % grab info for end time
    rad_idx = get_idx(x, s.rad(1));
    mou_idx = s.mou_idx(end-1);

    % levee complex
    fill([x_norm(1:rad_idx)/1000, fliplr(x_norm(1:rad_idx))/1000], ...
         [s.zed(1:rad_idx, end-1); baselinevect(1:rad_idx)], ...
         Cland1, 'EdgeColor', 'none'); % zed
    plot(x_norm(1:rad_idx)/1000, s.zed(1:rad_idx, end-1), ...
        'LineWidth', 2, 'Color', Cland2, 'LineStyle', '-'); % zed
    plot(x_norm(1:rad_idx)/1000, zed_i(1:rad_idx), 'LineWidth', 2, 'Color', Cland2, 'LineStyle', '--'); % initial 

    % lobe complex
    fill([x_norm(rad_idx:mou_idx+1)/1000, fliplr(x_norm(rad_idx:mou_idx+1))/1000], ...
         [s.zed(rad_idx:mou_idx+1, end); baselinevect(rad_idx:mou_idx+1)], ...
         Csubstr1, 'EdgeColor', 'none'); % zed
    plot(x_norm(rad_idx:mou_idx+1)/1000, s.zed(rad_idx:mou_idx+1, end), ...
        'LineWidth', 2, 'Color', Cland2, 'LineStyle', '-'); % zed

    % water in channel
    water = s.eta(:, end-1) + s.H(:, end-1);
    if true
        water(mou_idx:end) = 0;
    end
    fill([x_norm/1000, fliplr(x_norm)/1000], ...
         [water; baselinevect], ...
         Cwater1, 'EdgeColor', Cwater2, 'FaceAlpha', 0.6); % water
     
    % channel bed
    fill([x_norm/1000, fliplr(x_norm)/1000], ...
         [s.eta(:, end-1); baselinevect], ...
         Csubstr1,  'LineWidth', 1, 'LineStyle', '-', 'EdgeColor', Csubstr2); % channel bed line
    
    % original channel bed
    plot(x_norm/1000, eta_base, 'LineWidth', 2, 'Color', Csubstr2, 'LineStyle', '--')
     
    % loop bed colorlines for time
    for i = 1:length(pidx)
        p = pidx(i);
        plot(x_norm/1000, s.eta(:, p), 'LineWidth', 1, 'Color', colmap(i,:), 'LineStyle', '-'); % time lines
    end
     
    [t0i] = plot(x_norm/1000, eta_i, 'LineWidth', 1.5, 'Color', Csubstr2, 'LineStyle', '-'); % initial 
    [t0i] = plot(x_norm/1000, s.eta(:,end-1), 'LineWidth', 1, 'Color', Csubstr2, 'LineStyle', '-'); % final 
       
    % mouth and radius lines
    if s.rad(end-1) > s.rad(1)
        plot([x_norm(rad_idx_i)/1000, x_norm(rad_idx_i)/1000], [s.zed(rad_idx_i, end-1) s.eta(rad_idx_i, end-1)], ...
            'LineWidth', 1, 'Color', Cdgray, 'LineStyle', '-'); % initial radius line
    end
    plot([x_norm(rad_idx)/1000, x_norm(rad_idx)/1000], [s.zed(rad_idx, p) s.eta(rad_idx, p)], ...
        'LineWidth', 1, 'Color', 'k', 'LineStyle', '-'); % radius line
    plot([((mou_idx*dx)-s.rad(1))/1000, ((mou_idx*dx)-s.rad(1))/1000], ...
         [s.zed(mou_idx, p) s.eta(mou_idx, p)], ...
        'LineWidth', 1, 'Color', Csubstr2, 'LineStyle', '-');
    
    art_avul_idx = get_idx(x, (s.mou_idx(end-1)*dx - (data_struct.avulloc*d.Lblow)));
    text(0.6, 0.85, ['${L_D}^* = ', num2str(data_struct.avulloc), '$'], 'Units', 'Normalized')
    plot(x_norm(avul_pt)/1000, s.eta(avul_pt,end-1), ...
            'ko', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
        
    text(0.1, 0.85, char('a'+axidx-1), 'units', 'normalized', 'interpreter', 'latex', 'fontsize', 12)
        
    % clean up
    xlim([(d.gammapex-s.rad(1))/1000, max_x/1000]);
    ylim([-10 8])
    xlabel('dist. from delta coastline (km)','FontSize', 14);
    ylabel('elevation (m)','FontSize', 10);
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
    box on
    set(gca, 'Layer', 'top')
    hold off

end


function [pfig] = plot_avulsion_cycle_records(data_struct, pfig, axidx)

    s = data_struct.s;
    d = s.d;
    global CP

    x = d.x;
    dx = d.dx;
    Hnbf = d.Hnbf;
    L = max(d.x);
    nx = size(x, 2);
    
    % params
    upperlim = size(s.eta, 2)-1;
    nt = upperlim ./ 365;
    nprint = 10;
    
    % setup
    pidx = [1:floor(upperlim/nprint):upperlim, upperlim];
    colmap = parula(length(pidx));
    
    % colors
    Cwater1 = [195, 225, 255] ./ 255;
    Cwater2 = [35, 174, 216] ./ 255;
    Cland1 = [147, 236, 140] ./ 255;
    Cland2 = [123, 208, 117] ./ 255; 
    Csubstr1 = [196, 181, 132] ./ 255;
    Csubstr2 = [100, 63, 15] ./ 255;
    Cdgray = [0.5 0.5 0.5];
    
    % initial values
    eta_i = s.eta(:, 1);
    zed_i = s.zed(:, 1);
    cfg = yellowriver_configuration();
    [eta_base, ~, ~] = setup_initial_channel('new', floor((nx-1)/2), Hnbf, cfg.S0, nx-1, dx, ...
                                             repmat(cfg.Bc0, 1, nx), cfg.H0, cfg.Cf, ...
                                             cfg.Qwbf, cfg.eta_basin);
    
    % datums
    baselinevect = -20 * ones(nx, 1);
    gammapexxy = [cosd(d.gamm/2)*(d.gammapex), sind(d.gamm/2)*(d.gammapex)];
    
    x_norm = x - s.rad(1); % normalization x vector for rad_i for long profiles
    max_x = ((max(s.mou_idx)*dx)- s.rad(1))*1.5; % max x limit for long profiles
    
    [plan_i] = make_planform(x, d.gamm, d.gammapex, s.rad(1), false);
    xy_norm = [plan_i(end, 1), plan_i(1, 2)]; % normalization x and y scalar for planform plots
    [plan_i_norm] = plan_i - xy_norm;
    
    %% calculations
    % calculate the aggradation metrics
    agg = s.eta(:, 2:end) - s.eta(:, 1:end-1);
    rad_idx_i = get_idx(x, s.rad(1));
    avul_pt = s.avul(end);
    
    %% long prof
    sp1 = subplot(2, 3, axidx); % long profile
    hold on

    % grab info for end time
    rad_idx = get_idx(x, s.rad(1));
    mou_idx = s.mou_idx(end-1);

    % levee complex
    fill([x_norm(1:rad_idx)/1000, fliplr(x_norm(1:rad_idx))/1000], ...
         [s.zed(1:rad_idx, end-1); baselinevect(1:rad_idx)], ...
         Cland1, 'EdgeColor', 'none'); % zed
    plot(x_norm(1:rad_idx)/1000, s.zed(1:rad_idx, end-1), ...
        'LineWidth', 2, 'Color', Cland2, 'LineStyle', '-'); % zed
    plot(x_norm(1:rad_idx)/1000, zed_i(1:rad_idx), 'LineWidth', 2, 'Color', Cland2, 'LineStyle', '--'); % initial 

    % lobe complex
    fill([x_norm(rad_idx:mou_idx+1)/1000, fliplr(x_norm(rad_idx:mou_idx+1))/1000], ...
         [s.zed(rad_idx:mou_idx+1, end); baselinevect(rad_idx:mou_idx+1)], ...
         Csubstr1, 'EdgeColor', 'none'); % zed
    plot(x_norm(rad_idx:mou_idx+1)/1000, s.zed(rad_idx:mou_idx+1, end), ...
        'LineWidth', 2, 'Color', Cland2, 'LineStyle', '-'); % zed

    % water in channel
    water = s.eta(:, end-1) + s.H(:, end-1);
    if true
        water(mou_idx:end) = 0;
    end
    fill([x_norm/1000, fliplr(x_norm)/1000], ...
         [water; baselinevect], ...
         Cwater1, 'EdgeColor', Cwater2, 'FaceAlpha', 0.6); % water
     
    % channel bed
    fill([x_norm/1000, fliplr(x_norm)/1000], ...
         [s.eta(:, end-1); baselinevect], ...
         Csubstr1,  'LineWidth', 1, 'LineStyle', '-', 'EdgeColor', Csubstr2); % channel bed line
     
    % original channel bed
    plot(x_norm/1000, eta_base, 'LineWidth', 2, 'Color', Csubstr2, 'LineStyle', '--')
    
    % loop bed colorlines for time
    for i = 1:length(pidx)
        
        p = pidx(i);
        plot(x_norm/1000, s.eta(:, p), 'LineWidth', 1, 'Color', colmap(i,:), 'LineStyle', '-'); % time lines
    end
     
    [t0i] = plot(x_norm/1000, eta_i, 'LineWidth', 1.5, 'Color', Csubstr2, 'LineStyle', '-'); % initial 
    [t0i] = plot(x_norm/1000, s.eta(:,end-1), 'LineWidth', 1, 'Color', Csubstr2, 'LineStyle', '-'); % final 
       
    % mouth and radius lines
    if s.rad(end-1) > s.rad(1)
        plot([x_norm(rad_idx_i)/1000, x_norm(rad_idx_i)/1000], [s.zed(rad_idx_i, end-1) s.eta(rad_idx_i, end-1)], ...
            'LineWidth', 1, 'Color', Cdgray, 'LineStyle', '-'); % initial radius line
    end
    plot([x_norm(rad_idx)/1000, x_norm(rad_idx)/1000], [s.zed(rad_idx, p) s.eta(rad_idx, p)], ...
        'LineWidth', 1, 'Color', 'k', 'LineStyle', '-'); % radius line
    plot([((mou_idx*dx)-s.rad(1))/1000, ((mou_idx*dx)-s.rad(1))/1000], ...
         [s.zed(mou_idx, p) s.eta(mou_idx, p)], ...
        'LineWidth', 1, 'Color', 'k', 'LineStyle', '--');
    
    art_avul_idx = get_idx(x, (s.mou_idx(end-1)*dx - (data_struct.avulloc*d.Lblow)));
    plot(x_norm(art_avul_idx)/1000, s.eta(art_avul_idx,1), 'ko')
    text(0.5, 0.8, ['${L_D}^* = ', num2str(data_struct.avulloc), '$'], 'Units', 'Normalized')
    plot(x_norm(avul_pt)/1000, s.eta(avul_pt,end-1), ...
            'ko', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
        
    text(0.1, 0.85, char('a'+axidx-1), 'units', 'normalized', 'interpreter', 'latex', 'fontsize', 12)
        
    % clean up
    xlim([(d.gammapex-s.rad(1))/1000, max_x/1000]);
    ylim([-10 8])
    if axidx == 1
        ylabel('elevation (m)','FontSize', 10);
    end
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
    box on
    set(gca, 'Layer', 'top')
    hold off

    
    %% aggradation surface
    sp2 = subplot(2, 3, axidx+3); hold on;
    xx = 5; yy=3;
    [aggX, aggY] = meshgrid(x_norm./1000, (1:size(agg, 2))/365);
    [aggX] = BlockMean(aggX, xx, yy);
    [aggY] = BlockMean(aggY, xx, yy);
    [agg_ds] = BlockMean(agg', xx, yy)';
    surf(aggX, aggY, agg_ds'.*1000, 'EdgeColor', 'none');
    
    nt = size(agg,2);
    nyr = ceil(nt/365);
    plot3([0 0], [0 nt/365], [100 100], 'k-', 'LineWidth', 1)
    plot3((repmat(s.mou_idx(end-1), nyr+1,1)*dx-s.rad(1))/1000, [0:(nyr)]', repmat(100,nyr+1,1), ...
            'k--', 'LineWidth', 1) % mouth at end of TA'
    plot3(x_norm(avul_pt)/1000, nt/365, 100, ...
            'ko', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
    plot3(x_norm(art_avul_idx)/1000, 0, 100, 'ko')
    pause(0.5)
    
    text(0.1, 0.85, char('a'+axidx+3-1), 'units', 'normalized', 'interpreter', 'latex', 'fontsize', 12)
    
    xlabel('dist. from delta coastline (km)')
    ylabel('elapsed time (yr)')
    view([0 90])
    ylim([0 size(s.eta,2)/365])
    xlim([(d.gammapex-s.rad(1))/1000, max_x/1000]);
    colormap(flipud(redblue(21)))
    caxis([-1 1])
    
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
    set(gca, 'Layer', 'top')
    box on
    drawnow
    
    if axidx > 2
        orig_dims = get(sp2, 'position');
        cb = colorbar;
        ylabel(cb, 'deposition / erosion (mm/yr)')
        set(sp2, 'position', orig_dims);
        
    end
    
end