function movie_colored(movFig, data_path, nprint, save_path)

    if isstr(data_path)
        load(data_path)
    else
        warning('assuming data stucture of s and d instead of path')
        s = data_path.s;
        d = s.d;
        %cfg = data_path.cfg;
    end

    % convenience
    x = d.x;
    dx = d.dx;
    L = max(d.x);
    nx = size(x, 2);
    
    % params
    upperlim = size(s.eta, 2)-1;
    nt = upperlim ./ 365;
    
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
    Clgray = [0.8 0.8 0.8];
    Cdgray = [0.5 0.5 0.5];
    
    % initial values
    rad_idx_i = get_idx(x, s.rad(1));
    eta_i = s.eta(:, 1);
    zed_i = s.zed(:, 1);
    rad_i = s.rad(1)/1000;
    
    % datums
    baselinevect = -20 * ones(nx, 1);
    % timelinevect = ; % length needs to be set each time
    ymax = max(s.mou_idx)*dx/1000 * 1.1;
    ymin = nanmin(s.avul)*dx/1000 * 0.975;
    gammapexxy = [cosd(d.gamm/2)*(d.gammapex), sind(d.gamm/2)*(d.gammapex)];
%     rad_ixy = [cosd(d.gamm/2)*(s.rad(1)), sind(d.gamm/2)*(s.rad(1))];
    
    x_norm = x - s.rad(1); % normalization x vector for rad_i for long profiles
    max_x = max([((max(s.mou_idx)*dx)- s.rad(1))*1.5, 20e3]); % max x limit for long profiles
    
    [plan_i] = make_planform(x, d.gamm, d.gammapex, s.rad(1), false);
    xy_norm = [plan_i(end, 1), plan_i(1, 2)]; % normalization x and y scalar for planform plots
    [plan_i_norm] = plan_i - xy_norm;
    xy_baseline = min(plan_i_norm);
    gammapexxy_norm = gammapexxy - xy_norm;
    max_xy = [cosd(d.gamm/2)*(max(s.mou_idx)*dx), sind(d.gamm/2)*(max(s.mou_idx)*dx)]*1.1 - xy_norm; % max xy limit for planform profiles
    
    %%% FOR MANUAL X LIMS
%     max_x = 40e3;
%     max_xy = [cosd(d.gamm/2)*(s.rad(1)+max_x), sind(d.gamm/2)*(s.rad(1)+max_x)]*1.1 - xy_norm; % max xy limit for planform profiles
    
    figure(movFig);
    for f = 1:length(pidx)
        
        %% long prof
            subplot(2, 2, [1 2]); % long profile
            cla
            hold on
            
            p = pidx(f);
            
            % grad info for current line
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
                 Csubstr1, 'EdgeColor', 'none'); % zed
            plot(x_norm(rad_idx:mou_idx+1)/1000, s.zed(rad_idx:mou_idx+1, p), ...
                'LineWidth', 2, 'Color', Cland2, 'LineStyle', '-'); % zed
            
            % water in channel
            water = s.eta(:, p) + s.H(:, p);
            if true
                water(mou_idx:end) = 0;
            end
                
            fill([x_norm/1000, fliplr(x_norm)/1000], ...
                 [water; baselinevect], ...
                 Cwater1, 'EdgeColor', Cwater2, 'FaceAlpha', 0.6); % water
            
            % channel bed
            fill([x_norm/1000, fliplr(x_norm)/1000], ...
                 [s.eta(:, p); baselinevect], ...
                 Csubstr1,  'LineWidth', 1, 'LineStyle', '-', 'EdgeColor', Csubstr2); % channel bed line
            [t0f] = plot(x_norm/1000, eta_i, 'LineWidth', 2, 'Color', Csubstr2, 'LineStyle', '-'); % initial 
            
            % mouth and radius lines
            if s.rad(p) > s.rad(1)
                plot([x_norm(rad_idx_i)/1000, x_norm(rad_idx_i)/1000], [s.zed(rad_idx_i, p) s.eta(rad_idx_i, p)], ...
                    'LineWidth', 1, 'Color', Cdgray, 'LineStyle', '-'); % initial radius line
            end
            plot([x_norm(rad_idx)/1000, x_norm(rad_idx)/1000], [s.zed(rad_idx, p) s.eta(rad_idx, p)], ...
                'LineWidth', 1, 'Color', 'k', 'LineStyle', '-'); % radius line
            plot([((mou_idx*dx)-s.rad(1))/1000, ((mou_idx*dx)-s.rad(1))/1000], ...
                 [s.zed(mou_idx, p) s.eta(mou_idx, p)], ...
                'LineWidth', 1, 'Color', 'k', 'LineStyle', '--');
            
            % clean up
            xlim([(d.gammapex-s.rad(1))/1000, max_x/1000]);
            ylim([-10 8])
            xlabel('x distance (km)','FontSize', 14);
            ylabel('z distance (m)','FontSize', 14);
            set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
            box on
            set(gca, 'Layer', 'top')
        hold off
        
        %% planform
        subplot(2, 2, 3); % planform
            cla
            hold on

            % fill canvas with water
            fill([xy_baseline(1)/1000, xy_baseline(1)/1000, L/1000, L/1000], ...
                 [xy_baseline(2)/1000, L/1000, L/1000, xy_baseline(2)/1000], ...
                 Cwater1, 'EdgeColor', Cwater1, 'FaceAlpha', 0.6); % water

            % fill delta material
            [plan] = make_planform(x, d.gamm, d.gammapex, s.rad(p), true);
            plan_norm = plan - (xy_norm);
            fill(plan_norm(:, 1)/1000, plan_norm(:, 2)/1000, Cland1, 'EdgeColor', 'none')
            fill([xy_baseline(1)/1000, plan_norm(1, 1)/1000, plan_norm(1, 1)/1000, xy_baseline(1)/1000], ...
                 [plan_norm(1, 2)/1000, plan_norm(1, 2)/1000, xy_baseline(2)/1000, xy_baseline(2)/1000], Cland1, 'EdgeColor', 'none')
            fill([xy_baseline(1)/1000, plan_norm(end-1, 1)/1000, plan_norm(end-1, 1)/1000, xy_baseline(1)/1000], ...
                 [plan_norm(end-1, 2)/1000, plan_norm(end-1, 2)/1000, xy_baseline(2)/1000, xy_baseline(2)/1000], Cland1, 'EdgeColor', 'none')
            
            % delta edge line
            if s.rad(p) > s.rad(1)
                plot([xy_baseline(1)/1000; plan_i_norm(1:end, 1)/1000; plan_i_norm(end, 1)/1000], ...
                     [plan_i_norm(1, 2)/1000; plan_i_norm(1:end, 2)/1000; xy_baseline(2)/1000], ...
                     'Color', Cdgray)
            end
            plot([xy_baseline(1)/1000; plan_norm(1:end-1, 1)/1000; plan_norm(end-1, 1)/1000], ...
                 [plan_norm(1, 2)/1000; plan_norm(1:end-1, 2)/1000; xy_baseline(2)/1000], ...
                 'Color', 'k')
           
            % plot lobe
            p = pidx(f);
            [river, lobe] = make_lobeplan(x, d.gamm, s.rad(p), s.mou_idx(p), d.Bo0, dx);
            river_norm = river - xy_norm;
            lobe_norm = lobe - xy_norm;
            if s.mou_idx(p)*dx > s.rad(p)
                fill(lobe_norm(:,1)/1000, lobe_norm(:,2)/1000, Csubstr1, 'EdgeColor', 'none');
                plot(lobe_norm(:,1)/1000, lobe_norm(:,2)/1000, 'Color', 'k', 'LineStyle', '--');
            end
            
            % plot floodplain area
            fill([0, (cosd(d.gamm/2)*d.gammapex)/1000, fliplr([d.Bf0/2/1000, (d.Bf0/2 + cosd(d.gamm/2)*d.gammapex)/1000])]'-xy_norm/1000, ...
                 [d.Bf0/2/1000, (d.Bf0/2 + sind(d.gamm/2)*d.gammapex)/1000, fliplr([0 (sind(d.gamm/2)*d.gammapex)/1000])]'-xy_norm/1000, ...
                 Cland2, 'EdgeColor', 'none');

            plot(river_norm(:, 1)/1000, river_norm(:, 2)/1000, ...
                'LineWidth', 1.5, 'Color', Cwater2, 'LineStyle', '-'); % river centerline
            
            axis equal
%             xlim([0.3*L/1000 0.5*L/1000]);
%             ylim([0.3*L/1000 0.5*L/1000])
            xlim([gammapexxy_norm(1)/1000 max_xy(1)/1000]);
            ylim([gammapexxy_norm(2)/1000 max_xy(2)/1000])
            xlabel('x'' distance (km)','FontSize', 14);
            ylabel('y'' distance (km)','FontSize', 14);
            set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');

            box on
            set(gca, 'Layer', 'top')
        hold off
        
        %% data / discharge
        subplot(2, 2, 4)
        cla
        hold on
        if exist('cfg')
            if strcmp(cfg.Qw_str, 'historical')
                % Qing case
                qwdates = datenum('1976-01-01'):datenum('1996-12-31');
                newqw = [s.Qw, repmat(s.Qw(end), 1, 6)];
                plot(qwdates, newqw, 'Color', Cwater2, 'LineWidth', 1.25)
                plot(qwdates(pidx(f)), newqw(pidx(f)), 'k.', 'MarkerSize', 15)
                datetick('x')
                xlim([datenum('1975-12-31') datenum('1997-01-01')])
                ylabel('discharge (m^3/s)','FontSize', 14);
%                 set(gcf, 'Pos', [50 100 900 500], 'PaperPositionMode', 'auto')
            elseif strcmp(cfg.Qw_str, 'mean_historical')
                % longterm case
                nskip = 10;
                avul_idx = find(~isnan(s.avul(1:p)));
                timeseries_y = (1:nskip:p) ./ 365;
                timeseries_x = ymin-rad_i * ones(length(timeseries_y), 1)';
                fill([ymin-rad_i, ymax*1.1-rad_i, ymax*1.1-rad_i, ymin-rad_i], ...
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
                ylim([0, nt])
                xlim([ymin-rad_i, ymax-rad_i])
                ylabel('elapsed time (years)','FontSize', 14)
                xlabel('dist. from initial (km)','FontSize', 14);
                
            end
        end
            
            
            set(movFig, 'Pos', [50 100 650 650], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
            set(gca, 'Layer', 'top')
            box on
        
        %% update
        drawnow
        
        
        %% save it
        if ~isnan(save_path)   
            print('-dpng', '-r300', sprintf(['../movies/', save_path, '/%04d.png'], f));
            if p == 1
                print('-depsc', '-r300', '-painters', ['../figs/', save_path, '_frame1.eps']);
            end
            if p == upperlim
                print('-depsc', '-r300', '-painters', ['../figs/', save_path, '_framelast.eps']);
            end
        end
        
    end
    
   

end