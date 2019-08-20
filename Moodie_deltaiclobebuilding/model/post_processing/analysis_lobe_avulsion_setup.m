function analysis_lobe_avulsion_setup()

    clear all
    
    %% load data
    % add things to path
    source_path = genpath(fullfile('..', 'source'));
    output_data_path = genpath(fullfile('..', 'output_data'));
    plotting_path = genpath(fullfile('..', 'post_processing', 'plotting_scripts'));
    addpath(source_path, output_data_path, plotting_path);
    
    [CP] = color_palette();

    %% calculate a bunch of values for a single avulsion cycle
    if false
        load('longterm_data.mat')
        load('longterm_analysis.mat')
        
        avul_number = 14;
        time_idx = and(1:cnt.delap > analysis.avul_idx(avul_number), 1:cnt.delap <= analysis.avul_idx(avul_number+1));
        
        s = temporal_s_subset(s, time_idx);
        
        save('../output_data/longterm_data_nth.mat', 's', 'd', 'gammapex', 'Hnbf')

    else
        load('longterm_data_nth.mat')
        load('longterm_analysis.mat')
        
    end
    
    stat_summ = [(analysis.avul_time(15))./365, analysis.avul_len(15)./1000, analysis.lobe_len(15)./1000, analysis.avul_len(15)./analysis.Lblow, analysis.lobe_len(15)./analysis.Lblow];
    
    x = d.x;
    dx = d.dx;
    
    %% calculations
    % calculate the aggradation metrics
    agg = s.eta(:, 2:end) - s.eta(:, 1:end-1);
    agg_norm = bsxfun(@rdivide, agg, s.qs(1,2:end));
    rad_idx_i = get_idx(x, s.rad(1));
        
%     agg(rad_idx_i:end, :) = NaN;
%     agg(1:gammapex_idx, :) = NaN;
%     agg(agg<=0) = NaN;
%     agg(Qw > 1000) = NaN;
        
    [max_agg_val, max_agg_idx] = max(agg, [], 1);
    dist_upstream = (rad_idx_i - (max_agg_idx)) * dx / 1000;

    % aux variables
    lobe_length = (s.mou_idx(2:end)*dx - s.rad(2:end)) / 1000; % km

    avul_pt = s.avul(end);
    up_pt = avul_pt-30; % 20 km = 30 nodes
    down_pt = avul_pt+30;

    
    %% aggradation surface
    fig1 = figure(); hold on;
    
    [aggX, aggY] = meshgrid((x - s.mou_idx(end-1)*dx)./1000, (1:size(agg, 2))/365);
    agg_surf = surf(aggX, aggY, agg'.*1000, 'EdgeColor', 'none');
    
    xlabel('dist. from mouth at avul. time (km)')
    ylabel('elapsed time (yr)')
    view([0 90])
    ylim([0 size(agg,2)/365])
%     xlim([gammapex- s.mou_idx(end-1)*dx, ((max(s.mou_idx)*dx)- s.mou_idx(end-1)*dx)*1.5]./1000)
    xlim([-80 10])
    colormap(flipud(redblue(21)))
    cb = colorbar;
        ylabel(cb, 'deposition / erosion (mm/yr)')
    caxis([-1 1])
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
    set(gca, 'Layer', 'top')
    box on
    set(fig1, 'Pos', [100 100 400 400]);
    set(fig1, 'PaperPositionMode', 'auto')
    drawnow
    pause(1)
    set(gca,'Visible','off')
    print('-dpng', '-r300', '-opengl', '../figs/lobe_avulsion_aggerr.png');
    
    set(gca,'Visible','on')
    plot3([0 0], [0 size(agg,2)/365], [100 100], 'k-', 'LineWidth', 1.5)
    plot3((repmat(s.rad(1), ceil(size(agg,2)/365)+1,1)-s.mou_idx(end-1)*dx)/1000, ...
        [0:(ceil(size(agg,2)/365))]', ...
        repmat(100,ceil(size(agg,2)/365)+1,1), 'k--', 'LineWidth', 1.5)
    plot3((repmat(avul_pt,ceil(size(agg,2)/365)+1,1)*dx-s.mou_idx(end-1)*dx)/1000, ...
        [0:(ceil(size(agg,2)/365))]', ...
        repmat(100,ceil(size(agg,2)/365)+1,1), 'k:', 'LineWidth', 1.5)
    plot3((avul_pt*dx-s.mou_idx(end-1)*dx)/1000, ...
        size(agg,2)/365, ...
        100, 'k*', 'LineWidth', 1.5)
    drawnow
    print('-dpng', '-r300', '-opengl', '../figs/lobe_avulsion_aggerr_filled.png');
    delete(agg_surf)
    drawnow
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs', '../figs/lobe_avulsion_aggerr.eps'))
    
    
    
    %% aggradation profiles
    is_flood = s.Qw > 900;
    is_change = is_flood - [0 is_flood(1:end-1)];
    flood_starts = find(is_change == 1);
    flood_ends = find(is_change == -1);
    
    deta_floods = NaN(size(s.eta,1), length(flood_starts));
    deta_lows = NaN(size(s.eta,1), length(flood_starts));
    for i = 1:length(flood_starts)
        
        deta_floods(:,i) = s.eta(:, flood_ends(i)) - s.eta(:, flood_starts(i));
        
        if i >= 2
            deta_lows(:,i) = s.eta(:, +flood_starts(i)) - s.eta(:, flood_ends(i-1));
        end
        
    end
    
%     deta_lows(:,23) = NaN;

    deta_lows(rad_idx_i:end, :) = NaN;
    deta_floods(rad_idx_i:end, :) = NaN;
    
    x_norm = (x - s.mou_idx(end-1)*dx) ./ 1000;
%     x_limits = [gammapex- s.mou_idx(end-1)*dx, 0]./1000; % ((max(s.mou_idx)*dx)- s.rad(1))]./1000;
    x_limits = [-80 0];
    k=1;
    
    deta_lows_mean = nanmean(deta_lows(:,k:end), 2);
    deta_floods_mean = nanmean(deta_floods(:,k:end), 2);
    deta_lows_std = std(deta_lows(:,k:end), 0, 2, 'omitnan');
    deta_floods_std = std(deta_floods(:,k:end), 0, 2, 'omitnan');
    
    fig2 = figure(); hold on;
    plot(x_norm, deta_floods(:,k:end), '-', 'LineWidth', 1, 'Color', CP.c1_s2)
    plot(x_norm, deta_lows(:,k:end), '-', 'LineWidth', 1, 'Color', CP.c3_s2)
%     fill([x_norm, fliplr(x_norm)], ...
%          [deta_lows_mean+deta_lows_std; flipud(deta_lows_mean-deta_lows_std)]', ...
%          CP.c3_s2, 'EdgeColor', CP.c3_s3)
%     fill([x_norm, fliplr(x_norm)], ...
%          [deta_floods_mean+deta_floods_std; flipud(deta_floods_mean-deta_floods_std)]', ...
%          CP.c1_s2, 'EdgeColor', CP.c1_s3)
    plot(x_norm, deta_lows_mean, '-', 'LineWidth', 2, 'Color', CP.c3_s1)
    plot(x_norm, deta_floods_mean, '-', 'LineWidth', 2, 'Color', CP.c1_s1)
%     plot(x_norm, (s.eta(:, end-1) - s.eta(:, 1)), '-', 'LineWidth', 2, 'Color', CP.c4_s1)
    plot(x_norm, 1-((s.zed(:, end-1) - s.eta(:, end-1))./Hnbf), '-', 'LineWidth', 2, 'Color', CP.c4_s1)
    plot(x_norm, nanmean(deta_lows(:,k:end), 2)+nanmean(deta_floods(:,k:end), 2), '-', 'LineWidth', 2, 'Color', CP.c2_s1)
    
    plot([-100 20], [0 0], 'k--', 'LineWidth', 2)
    fill([s.rad(1)-s.mou_idx(end-1)*dx s.rad(1)-s.mou_idx(end-1)*dx 10 10]./1000, [-2 2 2 -2], CP.g3, 'EdgeColor', 'none')
    plot([s.rad(1)-s.mou_idx(end-1)*dx, s.rad(1)-s.mou_idx(end-1)*dx]./1000, [-1, 1], 'k-', 'LineWidth', 2)
%     plot([0, 0], [-0.5, 0.5], 'k:', 'LineWidth', 2)
    plot([((s.avul(end)*dx)-s.mou_idx(end-1)*dx), ((s.avul(end)*dx)-s.mou_idx(end-1)*dx)]./1000, [-1, 1], 'k:', 'LineWidth', 2)
    
    xlim(x_limits)
    ylim([-0.8, 0.8])
    xlabel('dist. from mouth at avul. time (km)')
    ylabel('normalized elevation (-)')
    
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
    set(gca, 'Layer', 'top')
    box on
    set(fig2, 'Pos', [100 100 400 250]);
    set(fig2, 'PaperPositionMode', 'auto')
    drawnow
    pause(1)
    print('-dpng', '-r300', '-opengl', '../figs/lobe_avulsion_setup.png');
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs', '../figs/lobe_avulsion_setup.eps'))
    
    
    
    %% do a version with splitting the data into chunks of time
    % resample to n= chunks of the run
    if true
        n=3;
        sect = floor(length(flood_starts)/n);
        deta_floods_chunk = NaN(size(s.eta, 1), n);
        deta_lows_chunk = NaN(size(s.eta, 1), n);
        for i = 1:n
            if i < n
                deta_floods_chunk(:,i) = nansum(deta_floods(:, ((i-1)*sect+1):((i-1)*sect + sect)), 2);
                deta_lows_chunk(:,i) = nansum(deta_lows(:, ((i-1)*sect+1):((i-1)*sect + sect)), 2);
            else
                deta_floods_chunk(:,i) = nansum(deta_floods(:, ((i-1)*sect+1):end), 2);
                deta_lows_chunk(:,i) = nansum(deta_lows(:, ((i-1)*sect+1):end), 2);
            end
        end
        deta_floods_chunk(rad_idx_i:end, :) = NaN;
        deta_lows_chunk(rad_idx_i:end, :) = NaN;
    end
    
    
    fig3 = figure(); hold on;
    % background lines
    se_thresh_line = plot([-100, 40], [0.5, 0.5], '-', 'LineWidth', 2, 'Color', CP.g3);
    avul_lcn_line = plot([((s.avul(end)*dx)-s.mou_idx(end-1)*dx), ((s.avul(end)*dx)-s.mou_idx(end-1)*dx)]./1000, [-1, 1], 'k:', 'LineWidth', 1.5);
    
    %plot each line with diff color for time (dark to light)
    f1_line = plot(x_norm, deta_floods_chunk(:,1)./Hnbf, '-', 'LineWidth', 2, 'Color', CP.c1_s1);
    l1_line = plot(x_norm, deta_lows_chunk(:,1)./Hnbf, '-', 'LineWidth', 2, 'Color', CP.c3_s1);
    f2_line = plot(x_norm, deta_floods_chunk(:,2)./Hnbf, '-', 'LineWidth', 2, 'Color', CP.c1_s3);
    l2_line = plot(x_norm, deta_lows_chunk(:,2)./Hnbf, '-', 'LineWidth', 2, 'Color', CP.c3_s3);
    f3_line = plot(x_norm, deta_floods_chunk(:,3)./Hnbf, '-', 'LineWidth', 2, 'Color', CP.c1_s2);
    l3_line = plot(x_norm, deta_lows_chunk(:,3)./Hnbf, '-', 'LineWidth', 2, 'Color', CP.c3_s2);

    se_line = plot(x_norm, 1-((s.zed(:, end-1) - s.eta(:, end-1))./Hnbf), '-', 'LineWidth', 2, 'Color', CP.c4_s1);
    agg_line = plot(x_norm, (sum(deta_floods_chunk,2)+sum(deta_lows_chunk,2))./Hnbf, '-', 'LineWidth', 2, 'Color', CP.c2_s1);
    
    plot([-100 20], [0 0], 'k--', 'LineWidth', 2)
    fill([s.rad(1)-s.mou_idx(end-1)*dx s.rad(1)-s.mou_idx(end-1)*dx 10 10]./1000, [-2 2 2 -2], CP.g3, 'EdgeColor', 'none')
    plot([s.rad(1)-s.mou_idx(end-1)*dx, s.rad(1)-s.mou_idx(end-1)*dx]./1000, [-1, 1], 'k-', 'LineWidth', 2)
    
    
    legend([f1_line l1_line f2_line l2_line f3_line l3_line se_line agg_line], ...
           {'1st third flood', '1st third low', '2nd third flood', '2nd third low', '3rd third flood', '3rd third low', ...
            'super-elevation', 'agg. during this avul. cycle'});
    xlabel('dist. from mouth at avul. time (km)')
    ylabel('normalized elevation (-)')
    xlim(x_limits)
    ylim([-0.8, 0.8])
    
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
    set(gca, 'Layer', 'top')
    box on
    set(fig3, 'Pos', [100 100 400 225]);
    set(fig3, 'PaperPositionMode', 'auto')
    drawnow
    pause(1)
    print('-dpng', '-r300', '-opengl', '../figs/lobe_avulsion_setup_chunks.png');
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs', '../figs/lobe_avulsion_setup_chunks.eps'))
    
    
    %% do a version with three "stacked" subplots of different lines compared
    fig5 = figure(); hold on;
    breaks = [0 1 2];
    c = 3;
    idxs = floor((1:c)*(length(flood_starts)-3)/c);
    
    plot(repmat([-100; 20], 1, length(breaks(2:end))), repmat(breaks(2:end), 2, 1)-0.5, '-', 'LineWidth', 2, 'Color', CP.g3) % divider lines (also setup threshold lines)
    for i = 1:length(breaks)
        i_idxs = (idxs(i)-2):(idxs(i)+2);
        deta_low_i_std = std(deta_lows(:, i_idxs), 0, 2, 'omitnan');
        deta_flood_i_std = std(deta_floods(:, i_idxs), 0, 2, 'omitnan');
        deta_low_i_maxmin = [max(deta_lows(:, i_idxs), [], 2); flipud(min(deta_lows(:, i_idxs), [], 2))];
        deta_flood_i_maxmin = [max(deta_floods(:, i_idxs), [], 2); flipud(min(deta_floods(:, i_idxs), [], 2))];
        
        [~, low_max_idx] = max(deta_lows(:, idxs(i)));
        [~, flood_min_idx] = min(deta_floods(:, idxs(i)));
        
        % plot standard deviations of +/-2 around the one
        nans_idx = isnan(deta_low_i_std);
        if true
            fill([x_norm(~nans_idx), fliplr(x_norm(~nans_idx))], ...
                 [breaks(i)+deta_floods(~nans_idx, idxs(i))+deta_flood_i_std(~nans_idx); flipud(breaks(i)+deta_floods(~nans_idx, idxs(i))-deta_flood_i_std(~nans_idx))]', ...
                 CP.c1_s2, 'EdgeColor', CP.c1_s3)
            fill([x_norm(~nans_idx), fliplr(x_norm(~nans_idx))], ...
                 [breaks(i)+deta_lows(~nans_idx, idxs(i))+deta_low_i_std(~nans_idx); flipud(breaks(i)+deta_lows(~nans_idx, idxs(i))-deta_low_i_std(~nans_idx))]', ...
                 CP.c3_s2, 'EdgeColor', CP.c3_s3)
        else
            fill([x_norm(~nans_idx), fliplr(x_norm(~nans_idx))], ...
                 breaks(i)+deta_flood_i_maxmin([~nans_idx; flipud(~nans_idx)]), ...
                 CP.c1_s2, 'EdgeColor', CP.c1_s3)
            fill([x_norm(~nans_idx), fliplr(x_norm(~nans_idx))], ...
                 breaks(i)+deta_low_i_maxmin([~nans_idx; flipud(~nans_idx)]), ...
                 CP.c3_s2, 'EdgeColor', CP.c3_s3)
        end
        
        % plot the main lines
        plot(x_norm, breaks(i)+deta_floods(:, idxs(i)), '-', 'LineWidth', 2, 'Color', CP.c1_s1)
        plot(x_norm, breaks(i)+deta_lows(:, idxs(i)), '-', 'LineWidth', 2, 'Color', CP.c3_s1)
        
        % plot the min and max points
        plot(x_norm(flood_min_idx), breaks(i)+deta_floods((flood_min_idx), idxs(i)), 'o', 'MarkerFaceColor', CP.c1_s1, 'MarkerEdgeColor', 'none')
        plot(x_norm(low_max_idx), breaks(i)+deta_lows(low_max_idx, idxs(i)), 'o', 'MarkerFaceColor', CP.c3_s1, 'MarkerEdgeColor', 'none')
        
        % plot the super elevation line
        plot(x_norm, breaks(i)+1-((s.zed(:, flood_starts(idxs(i))) - s.eta(:, flood_starts(idxs(i))))./Hnbf), '-', 'LineWidth', 2, 'Color', CP.c4_s1) % superelevation
        
        % plot the sum line
        plot(x_norm, breaks(i)+deta_floods(:, idxs(i))+deta_lows(:, idxs(i)), '-', 'LineWidth', 2, 'Color', CP.c2_s1) % sum of two lines
%         plot(x_norm, breaks(i)+ (s.eta(:, flood_ends(idxs(i)))-s.eta(:,1))./Hnbf, '-', 'LineWidth', 2, 'Color', CP.c2_s1) % sum of two lines
        
    end

    plot(repmat([-100; 20], 1, length(breaks)), repmat(breaks, 2, 1), 'k--', 'LineWidth', 2) % zeros lines
    fill([s.rad(1)-s.mou_idx(end-1)*dx s.rad(1)-s.mou_idx(end-1)*dx 10 10]./1000, [-5 5 5 -5], CP.g3, 'EdgeColor', 'none')
    plot([s.rad(1)-s.mou_idx(end-1)*dx, s.rad(1)-s.mou_idx(end-1)*dx]./1000, [-5, 5], 'k-', 'LineWidth', 1) % zone of prog line
    plot([((s.avul(end)*dx)-s.mou_idx(end-1)*dx), ((s.avul(end)*dx)-s.mou_idx(end-1)*dx)]./1000, [-5, 5], 'k:', 'LineWidth', 1.5) % avul loc line
    
    xlim(x_limits)
    ylim([breaks(1)-0.5, breaks(end)+0.5])
    xlabel('dist. from mouth at avul. time (km)')
    ylabel('normalized elevation (-)')
    
    orig_tick_locs = get(gca, 'YTick');
    orig_tick_vals = get(gca, 'YTickLabels');
    new_tick_vals = orig_tick_vals;
    new_tick_vals(2:2:end) = {'0'};
    new_tick_vals(3:2:end-1) = {'+/-0.5'};
    new_tick_vals(end) = {'+0.5'};
    set(gca, 'YTickLabels', new_tick_vals)
    
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
    set(gca, 'Layer', 'top')
    box on
    set(fig5, 'Pos', [100 100 400 400]);
    set(fig5, 'PaperPositionMode', 'auto')
    drawnow
    pause(1)
    print('-dpng', '-r300', '-opengl', '../figs/lobe_avulsion_setup_subs_stacked.png');
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs', '../figs/lobe_avulsion_setup_subs_stacked.eps'))
    
    %% do a version with three regular subplots of different lines compared
    fig6 = figure(); hold on;
    
    for i = 1:length(idxs)
        subplot(3, 1, length(idxs)+1-i); hold on;
            i_idxs = (idxs(i)-3):(idxs(i)+3);
            deta_low_i_std = std(deta_lows(:, i_idxs), 0, 2, 'omitnan');
            deta_flood_i_std = std(deta_floods(:, i_idxs), 0, 2, 'omitnan');
            deta_low_i_maxmin = [max(deta_lows(:, i_idxs), [], 2); flipud(min(deta_lows(:, i_idxs), [], 2))];
            deta_flood_i_maxmin = [max(deta_floods(:, i_idxs), [], 2); flipud(min(deta_floods(:, i_idxs), [], 2))];

            [~, low_max_idx] = max(deta_lows(:, idxs(i)));
            [~, flood_min_idx] = min(deta_floods(:, idxs(i)));

            % plot standard deviations of +/-3 around the one
            nans_idx = isnan(deta_low_i_std);
            if true
                fill([x_norm(~nans_idx), fliplr(x_norm(~nans_idx))], ...
                     [deta_floods(~nans_idx, idxs(i))+deta_flood_i_std(~nans_idx); flipud(deta_floods(~nans_idx, idxs(i))-deta_flood_i_std(~nans_idx))]', ...
                     CP.c1_s2, 'EdgeColor', CP.c1_s3)
                fill([x_norm(~nans_idx), fliplr(x_norm(~nans_idx))], ...
                     [deta_lows(~nans_idx, idxs(i))+deta_low_i_std(~nans_idx); flipud(deta_lows(~nans_idx, idxs(i))-deta_low_i_std(~nans_idx))]', ...
                     CP.c3_s2, 'EdgeColor', CP.c3_s3)
            else
                fill([x_norm(~nans_idx), fliplr(x_norm(~nans_idx))], ...
                     deta_flood_i_maxmin([~nans_idx; flipud(~nans_idx)]), ...
                     CP.c1_s2, 'EdgeColor', CP.c1_s3)
                fill([x_norm(~nans_idx), fliplr(x_norm(~nans_idx))], ...
                     deta_low_i_maxmin([~nans_idx; flipud(~nans_idx)]), ...
                     CP.c3_s2, 'EdgeColor', CP.c3_s3)
            end

            % plot the main lines
            plot(x_norm, deta_floods(:, idxs(i)), '-', 'LineWidth', 2, 'Color', CP.c1_s1)
            plot(x_norm, +deta_lows(:, idxs(i)), '-', 'LineWidth', 2, 'Color', CP.c3_s1)

            % plot the min and max points
            plot(x_norm(flood_min_idx), deta_floods((flood_min_idx), idxs(i)), 'o', 'MarkerFaceColor', CP.c1_s1, 'MarkerEdgeColor', 'none')
            plot(x_norm(low_max_idx), deta_lows(low_max_idx, idxs(i)), 'o', 'MarkerFaceColor', CP.c3_s1, 'MarkerEdgeColor', 'none')

            % plot the super elevation line
            plot(x_norm, 1-((s.zed(:, flood_starts(idxs(i))) - s.eta(:, flood_starts(idxs(i))))./Hnbf), '-', 'LineWidth', 2, 'Color', CP.c4_s1) % superelevation

            % plot the sum line
            plot(x_norm, deta_floods(:, idxs(i))+deta_lows(:, idxs(i)), '-', 'LineWidth', 2, 'Color', CP.c2_s1) % sum of two lines
    %         plot(x_norm, breaks(i)+ (s.eta(:, flood_ends(idxs(i)))-s.eta(:,1))./Hnbf, '-', 'LineWidth', 2, 'Color', CP.c2_s1) % sum of two lines

    
            plot([-100; 20], [0 0], 'k--', 'LineWidth', 2) % zeros line
            plot([-100; 20], [0.5 0.5], '-', 'LineWidth', 1.5, 'Color', CP.g4) % zeros line
            fill([s.rad(1)-s.mou_idx(end-1)*dx s.rad(1)-s.mou_idx(end-1)*dx 10 10]./1000, [-5 5 5 -5], CP.g4, 'EdgeColor', 'none')
            plot([s.rad(1)-s.mou_idx(end-1)*dx, s.rad(1)-s.mou_idx(end-1)*dx]./1000, [-5, 5], 'k-', 'LineWidth', 1) % zone of prog line
            plot([((s.avul(end)*dx)-s.mou_idx(end-1)*dx), ((s.avul(end)*dx)-s.mou_idx(end-1)*dx)]./1000, [-5, 5], 'k:', 'LineWidth', 1.5) % avul loc line
            
            plot([((s.mou_idx(flood_starts(idxs(i)))*dx-40e3)-s.mou_idx(end-1)*dx), ((s.mou_idx(flood_starts(idxs(i)))*dx-40e3)-s.mou_idx(end-1)*dx)]./1000, [-5, 5], '--', 'Color', CP.g2, 'LineWidth', 1.5) % backwater onset line

            xlim(x_limits)
            ylim([-0.6, 0.6])
            xlabel('dist. from mouth at avul. time (km)')
            ylabel('normalized elevation (-)')
            set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
            set(gca, 'Layer', 'top')
            box on
    end
    
    set(fig6, 'Pos', [100 100 400 530]);
    set(fig6, 'PaperPositionMode', 'auto')
    drawnow
    pause(1)
    print('-dpng', '-r300', '-opengl', '../figs/lobe_avulsion_setup_subs.png');
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs', '../figs/lobe_avulsion_setup_subs.eps'))
    
    
    
    
    %% make a plot of the time trajectories of avul location and nearby
    
    time_vect = 1:(size(s.eta,2)-1);
    avul_pt = avul_pt;
    up_pt = up_pt;
    down_pt = down_pt; % these get defined way at top
    
    se_avul = 1-((s.zed(avul_pt, 1:end-1) - s.eta(avul_pt, 1:end-1))./Hnbf);
    se_up = 1-((s.zed(up_pt, 1:end-1) - s.eta(up_pt, 1:end-1))./Hnbf);
    se_down = 1-((s.zed(down_pt, 1:end-1) -  s.eta(down_pt, 1:end-1))./Hnbf);
    
    ar_avul = (s.eta(avul_pt, 2:end) - s.eta(avul_pt, 1:end-1)) ./ 1; % aggradation rate
    ar_up = (s.eta(up_pt, 2:end) - s.eta(up_pt, 1:end-1)) ./ 1; 
    ar_down = (s.eta(down_pt, 2:end) - s.eta(down_pt, 1:end-1)) ./ 1; 
    
    fig4 = figure(); hold on;
    plot([-1, 30], [0.5, 0.5], 'k:', 'LineWidth', 2)
    avul_line = plot(time_vect./365, se_avul, '-', 'LineWidth', 2, 'Color', [0 0 0]);
    up_line = plot(time_vect./365, se_up, '-', 'LineWidth', 2, 'Color', CP.g2);
    down_line = plot(time_vect./365, se_down, '-', 'LineWidth', 2, 'Color', CP.g3);
    
    legend([avul_line, up_line, down_line], {'avul. location', '20 km upstream', '20 km downstream'}, 'Location', 'SouthEast')
    
    xlabel('elapsed time (yr)')
    ylabel('super-elevation (-)')
    xlim([0 (size(s.eta,2)-1)./365])
    ylim([0 0.6])
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
    set(gca, 'Layer', 'top')
    box on
    set(fig4, 'Pos', [100 100 400 175]);
    set(fig4, 'PaperPositionMode', 'auto')
    drawnow
    pause(1)
    print('-dpng', '-r300', '-opengl', '../figs/lobe_avulsion_setup_locs.png');
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs', '../figs/lobe_avulsion_setup_locs.eps'))
    
    
    fig5 = figure(); hold on;
    plot([-1, 30], [0.5, 0.5], 'k:', 'LineWidth', 2)
    avul_line = plot(time_vect./365, ar_avul, '-', 'LineWidth', 2, 'Color', [0 0 0]);
    up_line = plot(time_vect./365, ar_up, '-', 'LineWidth', 2, 'Color', CP.g2);
    down_line = plot(time_vect./365, ar_down, '-', 'LineWidth', 2, 'Color', CP.g3);
    
    legend([avul_line, up_line, down_line], {'avul. location', '20 km upstream', '20 km downstream'}, 'Location', 'NorthEast')
    
    xlabel('elapsed time (yr)')
    ylabel('aggradation rate (m/d)')
    xlim([0 6])
    ylim([-0.02 0.03])
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
    set(gca, 'Layer', 'top')
    box on
    set(fig5, 'Pos', [100 100 400 175]);
    set(fig5, 'PaperPositionMode', 'auto')
    drawnow
    pause(1)
    print('-dpng', '-r300', '-opengl', '../figs/lobe_avulsion_setup_locs_aggrate.png');
    print('-depsc', '-painters', '-r300', fullfile('..', 'figs', '../figs/lobe_avulsion_setup_locs_aggrate.eps'))
    
end

function calculate_se_for_all_avulsion_cycles()

    se = cell(length(avul_idxs)-1, 1);
    for i = 1:length(avul_idxs)-1
        time = [avul_idxs(i), avul_idxs(i+1)];
        Qw = s.Qw(time(1)+1:time(2));
        qs = s.qs((time(1)+1):time(2));
        
        % calculate the aggradation metrics
        agg = s.eta(:, (time(1)+1):time(2)) - s.eta(:, time(1):(time(2)-1));
        agg_norm = bsxfun(@rdivide, agg, qs);
        rad_idx_i = get_idx(x, s.rad(time(1)));
        
        agg(rad_idx_i:end, :) = NaN;
%         agg(1:gammapex_idx, :) = NaN;
%         agg(agg<=0) = NaN;
%         agg(Qw > 1000) = NaN;
        
        [max_agg_val, max_agg_idx] = max(agg, [], 1);
        dist_upstream = (rad_idx_i - (max_agg_idx)) * dx / 1000;
        
        % aux variables
        lobe_length = (s.mou_idx(time(1):time(2)-1)*dx - s.rad(time(1):time(2)-1)) / 1000; % km
        
        
        % do some averaging
        times = time(1):time(2)-1;
        [class] = categorical(s.mou_idx(times));
        cats = categories(class);
        avg_dist_upstream = zeros(length(cats),1);
        avg_lobe_length = zeros(length(cats),1);
        for j = 1:length(cats)
            idx = class == cats(j);
            avg_dist_upstream(j) = nanmean(dist_upstream(idx));
            avg_lobe_length(j) = str2num(cats{j}) - rad_idx_i;
        end
        
        se_i.max_agg_val = max_agg_val;
        se_i.max_agg_idx = max_agg_idx;
        se_i.dist_upstream = dist_upstream;
        se_i.avg_dist_upstream = avg_dist_upstream;
        se_i.avg_lobe_length = avg_lobe_length;
        se_i.lobe_length = lobe_length;
        se_i.Qw = Qw;
        se_i.rad_idx_i = rad_idx_i;
        se_i.time = time;
        se(i) = {se_i};
        
    end
    
    %% plot it
    cmap = parula(length(se));
    fig2 = figure();
    hold on;
    for i = 1:length(se)
        plot(se{i}.avg_lobe_length(1:end), se{i}.avg_dist_upstream(1:end), 'o', 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', 'k')
        leg{i} = [num2str(i)];
    end
    
end

