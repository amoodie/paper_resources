function analysis_longterm()

    clear all
    
    %% load data
    % add things to path
    source_path = genpath(fullfile('..', 'source'));
    output_data_path = genpath(fullfile('..', 'output_data'));
    plotting_path = genpath(fullfile('..', 'post_processing', 'plotting_scripts'));
    addpath(source_path, output_data_path, plotting_path);
    
    if false
        movFig = figure();
        movie_colored(movFig, 'longterm_data.mat', 600, 'longterm')
    end
    
    if false
        timeFig = figure();
        plot_longterm(timeFig, 'longterm_data.mat')
    end
    
%     load('longterm_data_noclosethresh_notrigger.mat')
%     load('longterm_analysis_noclosethresh_notrigger.mat')
    
%     load('longterm_data_05closethresh_notrigger.mat')
%     load('longterm_analysis_05closethresh_notrigger.mat')
    
%     load('longterm_data_noclosethresh_01trigger.mat')
%     load('longterm_analysis_noclosethresh_01trigger.mat')

%     load('longterm_data_newsetup.mat')
%     load('longterm_analysis_newsetup.mat')
    
    load('longterm_data.mat')
    load('longterm_analysis.mat')

    if false
        movie_longterm(s, d)
    end
    
    %% make series of regressions
    % an x axis to plot against in years
    timeseries_x = (1:size(s.eta, 2)) ./ 365;
    
    % determine the postspinup portion of the model
    navuls_spinup = 4; % number of avulsions to ignore
    radchange = [0, s.rad(2:end) - s.rad(1:end-1)];
    radchange_idx = find(radchange);
    postspinup_idx = (timeseries_x .* 365) > radchange_idx(navuls_spinup-1);
    summ_avuls_idx = (1:length(analysis.avul_time)) >= navuls_spinup;

    % make model of avulsion forward stepping
    avul_times = timeseries_x(and(~isnan(s.avul), postspinup_idx)); % years
    avul_locs = s.avul(and(~isnan(s.avul), postspinup_idx)) .* dx; % ms
    avul_mdl = fitlm(avul_times, avul_locs); % m/yr
    avul_rate = round(avul_mdl.Coefficients.Estimate(2)/1000, 2); % km/yr
    avul_rate_err = round(avul_mdl.Coefficients.SE(2)/1000, 3); % km/yr
    
    % bootstrap for error estimate of avulmdl
    for i = 1:30
        n = floor(length(avul_locs) / 3);
        X = datasample(avul_times, n, 'Replace', false)';
        Y = datasample(avul_locs, n, 'Replace', false)' ;
        btmdl = fitlm(X, Y);
        avul_rate_boot(i, 1) = btmdl.Coefficients.Estimate(2)/1000;
    end
    avul_rate_err_boot = std(avul_rate_boot);
    
    % small avulsions only rate
    small_avuls = analysis.avul_len(summ_avuls_idx) < 80e3;
    small_avul_mdl = fitlm(avul_times(small_avuls), avul_locs(small_avuls));
    small_avul_rate = small_avul_mdl.Coefficients.Estimate(2)/1000;
    
    
    % make model of delta prog rate
    delta_times = avul_times;
    delta_locs = s.rad(and(~isnan(s.avul), postspinup_idx));
    delta_mdl = fitlm(delta_times, delta_locs);
    delta_rate = round(delta_mdl.Coefficients.Estimate(2)/1000, 2);
    delta_rate_err = round(delta_mdl.Coefficients.SE(2)/1000, 3); % km/yr
    
    % bootstrap for error estimate on the radius progradation
    for i = 1:100
        n = floor(length(delta_locs) / 3);
        X = datasample(delta_times, n, 'Replace', false)';
        Y = datasample(delta_locs, n, 'Replace', false)' ;
        btmdl = fitlm(X, Y);
        delta_rate_boot(i, 1) = btmdl.Coefficients.Estimate(2)/1000;
    end
    delta_rate_err_boot = std(delta_rate_boot);
    % deltaradtable = array2table([timeseriesx(and(postspinup, radchange))', s.rad(and(postspinup, radchange))'], 'VariableNames', {'time', 'rad'});
    % deltaradtable.rad0 = deltaradtable.rad - deltaradtable.rad(1);
    % deltamdl = fitnlm(deltaradtable, 'rad0 ~ a * time ^(b)', [0 0]);
    % plot(timeseriesx(postspinup), ...
    %     (deltamdl.Coefficients.Estimate(1).* timeseriesx(postspinup) .^(deltamdl.Coefficients.Estimate(2)) + deltaradtable.rad(1))/1000, ...
    %     'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1.2)
    
    % plot bootstrapped values in histogram
    fig = figure();
        histogram(delta_rate_boot)
        ylabel(['count (n=', num2str(length(delta_rate_boot)), ')'])
        xlabel('bootstrapped rate (a=10, km/yr)')
%         ylim([0, 14])
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
        box on
        set(fig, 'Pos', [100 100 600 300]);
        set(fig, 'PaperPositionMode', 'auto')
        pause(2)
    
    % first half regression for delta prog
    hpt = floor(length(delta_locs)/2);
    btmdl_front = fitlm(delta_times(1:hpt), delta_locs(1:hpt));
    btmdl_back = fitlm(delta_times(hpt:end), delta_locs(hpt:end));
    btmdl_front.Coefficients.Estimate(2)/1000;
    btmdl_back.Coefficients.Estimate(2)/1000;
    

    %% plot the planform and rates figure
    fig2 = figure();
    subplot(1, 3, 1:2)
        cla; hold on;
        
        % plot the data 
        plot(timeseries_x, s.rad/1000, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 2)
        plot(timeseries_x, s.mou_idx*dx/1000, 'LineStyle', ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 2)
        plot(timeseries_x, s.avul*dx/1000, '*k')
        
%         plot(timeseries_x, s.Lb*dx/1000, 'LineStyle', ':', 'Color', [1 0.4 0.4], 'LineWidth', 2)
        

        % plot the avulsion regression
        plot(avul_times, (small_avul_mdl.Coefficients.Estimate(2) .* avul_times + small_avul_mdl.Coefficients.Estimate(1))/1000, ...
            'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1.2)
        text(0.8, 0.2, [num2str(round(small_avul_mdl.Coefficients.Estimate(2)/1000,3)) ' km/yr'], ...
            'Units', 'Normalized')
        
        % plot the delta progradation regression
        plot(avul_times, (delta_mdl.Coefficients.Estimate(2).* avul_times + delta_mdl.Coefficients.Estimate(1))/1000, ...
            'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1.2)
        text(0.2, 0.75, [num2str(round(delta_mdl.Coefficients.Estimate(2)/1000,3)) ' km/yr'], ...
            'Units', 'Normalized')
        
        xlabel('elapsed time (years)')
        ylabel('location from datum (km)')
%         legend('rad', 'mou', 'avulsion', 'Location', 'NorthWest')
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
        box on
        set(fig2, 'Pos', [100 100 600 300]);
        set(fig2, 'PaperPositionMode', 'auto')
        pause(2)
        print('-depsc', '-r300', '-painters', '../figs/longterm_planview.eps');
        print('-dpng', '-r200', '-opengl', '../figs/longterm_planview.png');

    
    %% plot the summary stats figure
    [CP] = color_palette();
    summ_avuls_idx = (1:length(analysis.avul_time)) >= navuls_spinup;
    fig3 = figure();
    subplot(1, 3, 1)
        cla; hold on;
            boxplot(analysis.avul_time(summ_avuls_idx)/365)
            plot(ones(size(analysis.avul_time(summ_avuls_idx))) .* 1.2 + randn(size(analysis.avul_time(summ_avuls_idx)))/75, ...
                analysis.avul_time(summ_avuls_idx)/365, 'ko', 'MarkerSize', 6)        
            ylabel('avulsion timescale (yr)')
    %         title('')
            set(gca, 'XTick', 1, 'XTickLabel', {'T_A'})
                    boxes = findobj(gca, 'tag', 'Box');
            xlim([0.8 1.3])
            ylim([min(analysis.avul_time(summ_avuls_idx)/365)*0.9, max(analysis.avul_time(summ_avuls_idx)/365)*1.1])
            medians = findobj(gca, 'tag', 'Median');
            lowerWhisker = findobj(gca, 'tag', 'Lower Whisker');
            upperWhisker = findobj(gca, 'tag', 'Upper Whisker');
            outliers = findobj(gca, 'tag', 'Outliers');
            lowerAdjVal = findobj(gca, 'tag', 'Lower Adjacent Value');
            upperAdjVal = findobj(gca, 'tag', 'Upper Adjacent Value');
            set(boxes, 'LineWidth', 2, 'Color', CP.c3_s3)
            set(medians, 'LineWidth', 2, 'Color', CP.c1_s1)
            set(lowerWhisker, 'LineWidth', 2)
            set(upperWhisker, 'LineWidth', 2)
            set(outliers, 'LineWidth', 2, 'Color', CP.c1_s1)
            set(lowerAdjVal, 'LineWidth', 2)
            set(upperAdjVal, 'LineWidth', 2)
            box on
            set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
            set(gca, 'Layer', 'top')
    subplot(1, 3, 2)
        cla; hold on;
            boxplot(analysis.avul_len(summ_avuls_idx)/1000)
            plot(ones(size(analysis.avul_len(summ_avuls_idx))).*1.2+randn(size(analysis.avul_time(summ_avuls_idx)))/75, ...
                analysis.avul_len(summ_avuls_idx)/1000, 'ko', 'MarkerSize', 6)
            ylabel('avulsion length (km)')
            set(gca, 'XTick', 1, 'XTickLabel', {'L_A'})
            xlim([0.8 1.3])
    %         ylim([26.4835 75])
            ylim([min(analysis.avul_len(summ_avuls_idx)/1000)*0.9, max(analysis.avul_len(summ_avuls_idx)/1000)*1.1])
            boxes = findobj(gca, 'tag', 'Box');
            medians = findobj(gca, 'tag', 'Median');
            lowerWhisker = findobj(gca, 'tag', 'Lower Whisker');
            upperWhisker = findobj(gca, 'tag', 'Upper Whisker');
            outliers = findobj(gca, 'tag', 'Outliers');
            lowerAdjVal = findobj(gca, 'tag', 'Lower Adjacent Value');
            upperAdjVal = findobj(gca, 'tag', 'Upper Adjacent Value');
            set(boxes, 'LineWidth', 2, 'Color', CP.c3_s3)
            set(medians, 'LineWidth', 2, 'Color', CP.c1_s1)
            set(lowerWhisker, 'LineWidth', 2)
            set(upperWhisker, 'LineWidth', 2)
            set(outliers, 'LineWidth', 2, 'Color', CP.c1_s1)
            set(lowerAdjVal, 'LineWidth', 2)
            set(upperAdjVal, 'LineWidth', 2)
            box on
            set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
            set(gca, 'Layer', 'top')
    subplot(1, 3, 3)
        cla; hold on;
            boxplot(analysis.lobe_len(summ_avuls_idx)/1000)
            plot(ones(size(analysis.lobe_len(summ_avuls_idx))).*1.2+randn(size(analysis.avul_time(summ_avuls_idx)))/75, ...
                analysis.lobe_len(summ_avuls_idx)/1000, 'ko', 'MarkerSize', 6)
            ylabel('lobe length (km)')
            set(gca, 'XTick', 1, 'XTickLabel', {'L_L'})
            xlim([0.8 1.3])
%             ylim([0 10])
            ylim([min(analysis.lobe_len(summ_avuls_idx)/1000)*0.9, max(analysis.lobe_len(summ_avuls_idx)/1000)*1.1])
            boxes = findobj(gca, 'tag', 'Box');
            medians = findobj(gca, 'tag', 'Median');
            lowerWhisker = findobj(gca, 'tag', 'Lower Whisker');
            upperWhisker = findobj(gca, 'tag', 'Upper Whisker');
            outliers = findobj(gca, 'tag', 'Outliers');
            lowerAdjVal = findobj(gca, 'tag', 'Lower Adjacent Value');
            upperAdjVal = findobj(gca, 'tag', 'Upper Adjacent Value');
            set(boxes, 'LineWidth', 2, 'Color', CP.c3_s3)
            set(medians, 'LineWidth', 2, 'Color', CP.c1_s1)
            set(lowerWhisker, 'LineWidth', 2)
            set(upperWhisker, 'LineWidth', 2)
            set(outliers, 'LineWidth', 2, 'Color', CP.c1_s1)
            set(lowerAdjVal, 'LineWidth', 2)
            set(upperAdjVal, 'LineWidth', 2)
            box on
            set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
            set(gca, 'Layer', 'top')
        set(fig3, 'Pos', [100 100 600 300]);
        set(fig3, 'PaperPositionMode', 'auto')
        print('-depsc', '-r300', '-painters', '../figs/longterm_record.eps');
        print('-dpng', '-r200', '-opengl', '../figs/longterm_record.png');
        
        
        %% calcualte coefficient of variation
        Qw1 = s.Qw(1:365);
        H1 = s.H(1, 1:365);
        
        Qw2 = s.Qw(365:(365*2));
        H2 = s.H(1, 365:(365*2));
        
        cv1 = std(H1) ./ mean(H1);
        cv2 = std(H2) ./ mean(H2);
end