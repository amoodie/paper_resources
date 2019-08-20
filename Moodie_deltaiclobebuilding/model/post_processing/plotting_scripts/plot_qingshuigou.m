function plot_qingshuigou(qingFig, s, d)
    
    %% setup params etc
    Tdays = size(s.eta, 2);

    % regression portion
    Y = s.mou_idx(s.mou_idx.*d.dx>s.rad(1)) .* d.dx;
    X = find(s.mou_idx.*d.dx>s.rad(1), 1, 'first'):Tdays;
    mdl = fitlm(X, Y);
    relcoords = [0.8, 0.3];
    CI = coefCI(mdl, 0.05);
    err = mean(abs(CI(2,:) - mdl.Coefficients.Estimate(2)));
    
    nprint = 10;
    upperlim = size(s.eta, 2);
    pidx = 1:floor(upperlim/nprint):upperlim;
    pidx = pidx(1:end-1);
    colmap = parula(length(pidx));

    end_rad_idx = find(d.x > s.rad(end), 1, 'first');
    
    p = Tdays;
    
    % convenience
    x = d.x;
    dx = d.dx;
    L = max(d.x);
    nx = size(x, 2);
    nt = length(s.rad) ./ 365;
    
    % colors
    Cwater1 = [195, 225, 255] ./ 255;
    Cwater2 = [35, 174, 216] ./ 255;
    Cland1 = [147, 236, 140] ./ 255;
    Cland2 = [123, 208, 117] ./ 255; 
    Csubstr1 = [196, 181, 132] ./ 255;
    Csubstr2 = [100, 63, 15] ./ 255;
    Clgray = [0.8 0.8 0.8];
    Cdgray = [0.5 0.5 0.5];
    
    % initial values
    rad_idx_i = get_idx(x, s.rad(1));
    [plan_i] = make_planform(x, d.gamm, d.gammapex, s.rad(1), false);
    eta_i = s.eta(:, 1);
    zed_i = s.zed(:, 1);
    rad_i = s.rad(1)/1000;
    
    % datums
    baselinevect = -20 * ones(nx, 1);
    xmax = max(s.mou_idx)*dx/1000 * 1.05;
    xmin = -10;
    gammapexxy = [cosd(d.gamm/2)*(d.gammapex), sind(d.gamm/2)*(d.gammapex)];
    
    x_norm = x - s.rad(1); % normalization x vector for rad_i for long profiles
    max_x = ((max(s.mou_idx)*dx)- s.rad(1))*1.5; % max x limit for long profiles
    
    [plan_i] = make_planform(x, d.gamm, d.gammapex, s.rad(1), false);
    xy_norm = [plan_i(end, 1), plan_i(1, 2)]; % normalization x and y scalar for planform plots
    [plan_i_norm] = plan_i - xy_norm;
    xy_baseline = min(plan_i_norm);
    gammapexxy_norm = gammapexxy - xy_norm;
    max_xy = [cosd(d.gamm/2)*(max(s.mou_idx)*dx), sind(d.gamm/2)*(max(s.mou_idx)*dx)]*1.1 - xy_norm; % max xy limit for planform profiles
    
    
    figure(qingFig);
    %% make the long profile view
    subplot(2, 2, [1 2]); % flux plot
    cla
    hold on
      
        % grab info for current line
        rad_idx = get_idx(x, s.rad(p));
        mou_idx = s.mou_idx(p);

        % levee complex
        fill([x_norm(1:rad_idx)/1000, fliplr(x_norm(1:rad_idx))/1000], ...
             [s.zed(1:rad_idx, p); baselinevect(1:rad_idx)], ...
             Cland1, 'EdgeColor', 'none'); % zed
        plot(x_norm(1:rad_idx)/1000, s.zed(1:rad_idx, p), ...
            'LineWidth', 2, 'Color', Cland2, 'LineStyle', '-'); % zed
        plot(x_norm(1:rad_idx)/1000, zed_i(1:rad_idx), 'LineWidth', 2, 'Color', Cland2, 'LineStyle', '--'); % initial 

        % lobe complex
        fill([x_norm(rad_idx:mou_idx+1)/1000, fliplr(x_norm(rad_idx:mou_idx+1))/1000], ...
             [s.zed(rad_idx:mou_idx+1, p); baselinevect(rad_idx:mou_idx+1)], ...
             Cland1, 'EdgeColor', 'none'); % zed
        plot(x_norm(rad_idx:mou_idx+1)/1000, s.zed(rad_idx:mou_idx+1, p), ...
            'LineWidth', 2, 'Color', Cland2, 'LineStyle', '-'); % zed

        % water in channel
        water_end = s.eta(:, p) + s.H(:, p);
        water_high = s.eta(:, 7550) + s.H(:, 7550);
        if true
            water_end(mou_idx:end) = 0;
            water_high(mou_idx:end) = 0;
        end
        fill([x_norm/1000, fliplr(x_norm)/1000], ...
             [water_end; baselinevect], ...
             Cwater1, 'EdgeColor', Cwater2, 'FaceAlpha', 0.6); % end water
        plot(x_norm/1000, water_high, 'LineWidth', 1.5, 'Color', [1,0.58,0.58]) % a late-stage high water flow

        % channel bed
        fill([x_norm/1000, fliplr(x_norm)/1000], ...
             [s.eta(:, p); baselinevect], ...
             Csubstr1,  'LineWidth', 1, 'LineStyle', '-', 'EdgeColor', 'none'); % channel bed line
        for pp = 1:length(pidx)
            % curr_mou_chg_idx = (s.eta(:, pp) - s.eta(:, 1)) > 0.01;
            ppp = pidx(pp);
            plot(x_norm/1000, s.eta(:, ppp), 'LineWidth', 1.5, 'Color', colmap(pp, :), 'LineStyle', '-'); % time lines
        end
        [t0f] = plot(x_norm/1000, eta_i, 'LineWidth', 2, 'Color', Csubstr2, 'LineStyle', '-'); % initial
        [t0f] = plot(x_norm/1000, s.eta(:, end), 'LineWidth', 1, 'Color', Csubstr2, 'LineStyle', '-'); % final
        

        % mouth and radius lines
        if s.rad(p) > s.rad(1)
            plot([x_norm(rad_idx_i)/1000, x_norm(rad_idx_i)/1000], [s.zed(rad_idx_i, p) s.eta(rad_idx_i, p)], ...
                'LineWidth', 1, 'Color', Cdgray, 'LineStyle', '-'); % initial radius line
        end
        plot([x_norm(rad_idx)/1000, x_norm(rad_idx)/1000], [s.zed(rad_idx, p) s.eta(rad_idx, p)], ...
            'LineWidth', 1, 'Color', 'k', 'LineStyle', '-'); % radius line
        plot([mou_idx*dx/1000, mou_idx*dx/1000], [s.zed(mou_idx, p) s.eta(mou_idx, p)], ...
            'LineWidth', 1, 'Color', 'k', 'LineStyle', '--');

        % clean up
        xlim([(d.gammapex-s.rad(1))/1000, max_x/1000]);
        ylim([-8 8])
        xlabel('x distance (km)','FontSize', 14);
        ylabel('z distance (m)','FontSize', 14);
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        box on
        set(gca, 'Layer', 'top')
        hold off
    
    %% make the planview with regression
    s4 = subplot(2, 2, 3);
    cla
    hold on
        
        nskip = 10;
        avul_idx = find(~isnan(s.avul(1:p)));
        timeseries_y = (1:nskip:p) ./ 365;
        timeseries_x = xmin-rad_i * ones(length(timeseries_y), 1)';
        fill([xmin-rad_i, xmax*1.1-rad_i, xmax*1.1-rad_i, xmin-rad_i], ...
              [0, 0, p/365, p/365], ...
              Cwater1, 'EdgeColor', 'none') % water area
        fill([s.mou_idx(1:nskip:p)*dx/1000-rad_i, timeseries_x], ... 
             [timeseries_y, fliplr(timeseries_y)], ...
             Csubstr1, 'EdgeColor', Csubstr2) % lobe area
        fill([s.rad(1:nskip:p)/1000-rad_i, timeseries_x], ...
             [timeseries_y, fliplr(timeseries_y)], ...
             Cland1, 'EdgeColor', Cland2) % delta inside
        fill([s.rad(1)-d.Lblow, s.mou_idx(avul_idx)*dx-d.Lblow, fliplr(s.mou_idx(avul_idx-1)*dx-d.Lblow) s.rad(1)]./1000-rad_i, ...
             [0 avul_idx, fliplr(avul_idx) 0]./365, ...
             [150 150 150]./255, 'EdgeColor', 'none', 'FaceAlpha', 0.4)
        plot(s.avul(1:p)*dx/1000-rad_i, (1:p)./365, '*k')
         
        % regression plot
        plot(((mdl.Coefficients.Estimate(2).*[find(s.mou_idx*d.dx>s.rad(1), 1, 'first'):Tdays] + mdl.Coefficients.Estimate(1))/1000)-rad_i, ...
             [find(s.mou_idx*d.dx>s.rad(1), 1, 'first'):Tdays]./365, ...
             'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1.2)
        
        % annualized plot
        plot([s.mou_idx(1), s.mou_idx(end)]*dx/1000 - rad_i, ...
             [0, Tdays]./365, ...
             'LineStyle', ':', 'Color', [0 0 0], 'LineWidth', 1.2)
        
        % annotation follows
        params = {['rate = ', num2str(round(mdl.Coefficients.Estimate(2)*365.25/1000, 2)), '\pm', num2str(round(err*365.25/1000, 2)), ' km/yr'], ...
            ['r^2 = ', num2str(round(mdl.Rsquared.Ordinary, 2))]};
        format = sprintf('%s\n', params{:});
        annot = text(relcoords(1), relcoords(2), format(1:end-1), ...
            'Color', [0 0 0], 'Parent', s4, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        
        params = {['rate = ', num2str(round((s.mou_idx(end)*dx/1000 - rad_i)/(Tdays/365), 2)), ' km/yr']};
        format = sprintf('%s\n', params{:});
        annot = text(0.1, 0.8, format(1:end-1), ...
            'Color', [0 0 0], 'Parent', s4, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        
        ylim([0, nt-1])
        xlim([xmin, xmax-rad_i])
        ylabel('elapsed time (years)','FontSize', 14)
        xlabel('dist. from initial (km)','FontSize', 14);
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(gca, 'Layer', 'top')
        box on
        
    set(gcf, 'Pos', [50 100 650 650], 'PaperPositionMode', 'auto')
    % print('-depsc', '-r300', '-painters', '../figs/qingshuigou.eps');
    % print('-dpng', '-r200', '-opengl', '../figs/qingshuigou.png');


    

end