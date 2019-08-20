function analysis_qingshuigou()

    % add things to path
    source_path = genpath(fullfile('..', 'source'));
    output_data_path = genpath(fullfile('..', 'output_data'));
    plotting_path = genpath(fullfile('..', 'post_processing', 'plotting_scripts'));
    addpath(source_path, output_data_path, plotting_path);
    
    if false
        movFig = figure();
        movie_qingshuigou(movFig, 'qingshuigou_data.mat')
    end
    
    if false
        movFig = figure();
        movie_colored(movFig, 'qingshuigou_data.mat',  600, 'qing')
    end
    
    load('qingshuigou_data.mat')
    
    s.Be = double(s.Be);
    s.Bc = double(s.Bc);
    
    % make the profile/planform qingshuigou plot
    if true
        qingFig = figure();
        plot_qingshuigou(qingFig, s, d)
    end


    % make the plot of mass conservation comparison through time
    if true
        [mc, calc_idx] = calculate_mass_conservation(s, Bc, Be, Bc0, Bf0, Bo0, rad_idx, sig, lastavul_idx, dx, con);
        plot_mass_conservation(mc, calc_idx)
    end
    s.v.Vol_inside_rad = mc.vol_in_rad_vect;
    s.v.Vol_inside_floodplain = mc.vol_floodplain_vect;
    s.v.Vol_lobe = mc.vol_to_lobe_vect;
    s.v.integrated_input = mc.vol_input_vect;

    [s] = volume_change_per_section(s, d, sig, Bc0, Bf0, Bo0, con);
    
    [s] = water_discharge_binned_plots(s, d);
    
%     Tdays = cnt.delap;
%     figure(); hold on
%     plot(1:Tdays, s.f.normcumtop , 'LineWidth', 1.2)
%     plot(1:Tdays, s.f.normcumfore , 'LineWidth', 1.2)
%     plot(1:Tdays, (s.f.lobecumkm3 ./ (s.f.topcumkm3 + s.f.lobecumkm3)), 'LineWidth', 1.2)
%     plot(1:Tdays, (s.f.topcumkm3 ./ (s.f.topcumkm3 + s.f.lobecumkm3)), 'LineWidth', 1.2)
%     legend('slope - top', 'slope - fore', 'fixed - top', 'fixed - lobe')

    
    %% determine the cause of the shutdown in progradation
    cum_volume_int = cumtrapz(s.qs(rad_idx, :) .* s.Bc(rad_idx, :))  * 86400;
    figure; 
    plot(cum_volume_int * 1e-9)
    title('cumulative flux of sediment past radius plane')
    xlabel('time (d)')
    ylabel('volume (km^3)')
    
    if false
        figure()
        dqsdx = s.qs(2:end, :) - s.qs(1:end-1, :);
    %     dqsdx(abs(dqsdx) < 0.0005) = NaN;
        surf(dqsdx, 'EdgeColor', 'None')
        view([0 90])
        caxis([-0.000005, 0.000005])
        title('sediment flux map')
        colorbar
        xlabel('time (d)')
        ylabel('idx, dist downstream *dx')
    end
    
    %% plot of superelevation
     % barplot of where deposition
%     tickidxs = 1:5:length(q.Qcnt);
    max_x = ((mou_idx*dx)- rad)*1.5; % max x limit for long profiles
    superelev = 1 - ( (zed-eta)./Hnbf );
    superelev(mou_idx:end) = NaN;
    f = figure();
%         bar(1:length(q.Qcnt), [q.delta_vol.alltP, q.lobe_vol.alltP], 1)
        plot((x-rad)/1000, superelev, '-', 'Color', [204, 0, 102]./255, 'LineWidth', 1.2)
        ylim([0 1])
        xlim([(d.gammapex-rad)/1000, max_x/1000]);
        xlabel('x distance (km)')
        ylabel('super-elevation (H_{bf})')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(f, 'Position', [500 500 450 250])
        set(f, 'PaperPositionMode', 'auto')
        drawnow
        print('-depsc', '-painters', '-r300',  '../figs/qingshuigou_superelevation.eps');
    
end
    
function [s] = water_discharge_binned_plots(s, d)
    
    % make the counted bins for the following plots etc
    [q.Qcnt, q.Qedg, q.Qbin] = histcounts(s.Qw, 50);
    [q.Qedgelabs] = makelabels(q.Qedg);
    
    % preallocate the table to add rows to
    q.delta_vol = cell2table(cell(0, 3), 'VariableNames', {'allt', 'pre', 'post'});
    q.lobe_vol = cell2table(cell(0, 3), 'VariableNames', {'allt', 'pre', 'post'});
    q.tot_vol = cell2table(cell(0, 3), 'VariableNames', {'allt', 'pre', 'post'});
    
    % index of dates before
    preMouidx = s.mou_idx == s.mou_idx(1); 
    
    % index of locations inside the delta
    rad_idx = get_idx(d.x, s.rad(1));
    delta_idx = logical( 1:size(s.eta,1) < rad_idx );
    
    % use the deposition only Vdtdx from volume_change_per_section
    Vdtdx = s.f.Vdtdx_dep;
    
    % cumulate the proportions in each for each time slice
    for b = 1:length(q.Qcnt)
        Qidx = q.Qbin == b;
        
        % delta volume deposition
        q.delta_vol(b, :) = num2cell( [sum(sum(Vdtdx(delta_idx, Qidx))), ...
                                       sum(sum(Vdtdx(delta_idx, and(Qidx, preMouidx)))), ...
                                       sum(sum(Vdtdx(delta_idx, and(Qidx, ~preMouidx))))] );
        
        % lobe volume deposition
        q.lobe_vol(b, :) = num2cell( [sum(sum(Vdtdx(~delta_idx, Qidx))), ...
                                      sum(sum(Vdtdx(~delta_idx, and(Qidx, preMouidx)))), ...
                                      sum(sum(Vdtdx(~delta_idx, and(Qidx, ~preMouidx))))] );
        
        % total volume deposition is sum of two
        q.tot_vol(b, :) = num2cell( [plus(q.delta_vol.allt(b), q.lobe_vol.allt(b)), ...
                                     plus(q.delta_vol.pre(b), q.lobe_vol.pre(b)), ...
                                     plus(q.delta_vol.post(b), q.lobe_vol.post(b))] );
        
        % calculate qs upstream input per second for each discharge
        q.qs_inputs(b) = sum( s.qs(1, Qidx) .* s.Bc(1, Qidx) .* 86400 );
        q.qs_radflux(b) = sum( s.qs(rad_idx, Qidx) .* s.Bc(rad_idx, Qidx) .* 86400 );

%         % calculate input qs per second for each discharge
%         idx = 1; % input
%         q.qs_inputs(b) = integrate_qs_flux(s.qs(:, Qidx), s.Bc(:, Qidx), idx);
%         
%         % calculate to lobe flux qs per second for each discharge
%         idx = rad_idx; % radius
%         q.qs_inputs(b) = integrate_qs_flux(s.qs(:, Qidx), s.Bc(:, Qidx), idx);
% 
%         % calculate input qs per second for each discharge
%         idx = 1; % input
%         q.qs_inputs(b) = integrate_qs_flux(s.qs(:, Qidx), s.Bc(:, Qidx), idx);

    end
    
    % append proportional columns to end
    q.delta_vol = [q.delta_vol, array2table([q.delta_vol{:,:} ./ q.tot_vol{:,:}], 'VariableNames', {'alltP', 'preP', 'postP'})];
    q.lobe_vol = [q.lobe_vol, array2table([q.lobe_vol{:,:} ./ q.tot_vol{:,:}], 'VariableNames', {'alltP', 'preP', 'postP'})];
    
    
    %% plots
    % bins of discharge
    f = figure();
        histogram(s.Qw, 0:1000:6000, 'Normalization', 'probability')
        vals = histcounts(s.Qw, 0:1000:6000, 'Normalization', 'probability');
        text(200:1000:6000, vals+0.1, arrayfun(@(x) num2str(round(x,3)), vals, 'Unif', 0))
        ylim([0 1])
        xlabel('water discharge (m^3/s)')
        ylabel('fraction')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(f, 'Position', [300 200 450 250])
        set(f, 'PaperPositionMode', 'auto')
%         print('-dpng', '-r300', '../figs/qingshuigou_discharge_binned.png');
        drawnow; pause(1)
        print('-depsc', '-painters', '-r300',  '../figs/qingshuigou_discharge_binned.eps');
        
    % bins of sediment discharge
    total_sed = sum(s.qs(1,:));
    sed_bins = 0:1000:6000;
    prop_sed = zeros(length(sed_bins)-1, 1);
    for i = 2:length(sed_bins)
        qidx = and( (s.Qw < sed_bins(i)), (s.Qw > sed_bins(i-1)) );
        prop_sed(i-1) = sum(s.qs(1, qidx)) ./ total_sed;
    end
    f = figure();
        bar(500:1000:5500, prop_sed, 1)
        xticks(0:1000:6000)
        text(200:1000:6000, prop_sed+0.1, arrayfun(@(x) num2str(round(x,3)), prop_sed, 'Unif', 0))
        ylim([0 1])
        xlim([-250, 6250])
        xlabel('water discharge (m^3/s)')
        ylabel('fraction of sediment')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(f, 'Position', [300 200 450 250])
        set(f, 'PaperPositionMode', 'auto')
%         print('-dpng', '-r300', '../figs/qingshuigou_sedimentdischarge_binned.png');
        drawnow
        pause(1)
        print('-depsc', '-painters', '-r300',  '../figs/qingshuigou_sedimentdischarge_binned.eps');
    
    % barplot of where deposition
    tickidxs = 1:5:length(q.Qcnt);
    f = figure();
        bar(1:length(q.Qcnt), [q.delta_vol.alltP, q.lobe_vol.alltP], 1)
        ylim([0 1])
        set(gca, 'XTick', tickidxs)
        set(gca, 'xTickLabels', q.Qedgelabs(tickidxs), 'xTickLabelRotation', 45)
        title('allt')
        legend('delta', 'lobe')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(f, 'Position', [500 500 450 250])
        set(f, 'PaperPositionMode', 'auto')
    
    % barplot of where during: all, pre, and post
    figure();
    subplot(3, 1, 1)
        bar(1:length(q.Qcnt), [q.delta_vol.alltP, q.lobe_vol.alltP])
        set(gca, 'XTick', tickidxs)
        set(gca, 'xTickLabels', q.Qedgelabs(tickidxs), 'xTickLabelRotation', 45)
        title('allt')
        legend('topset', 'foreset')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
    subplot(3, 1, 2)
        bar(1:length(q.Qcnt), [q.delta_vol.preP, q.lobe_vol.preP])
        set(gca, 'XTick', tickidxs)
        set(gca, 'xTickLabels', q.Qedgelabs(tickidxs), 'xTickLabelRotation', 45)
        title('pre')
        legend('topset', 'foreset')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
    subplot(3, 1, 3)
        bar(1:length(q.Qcnt), [q.delta_vol.postP, q.lobe_vol.postP])
        set(gca, 'XTick', tickidxs)
        set(gca, 'xTickLabels', q.Qedgelabs(tickidxs), 'xTickLabelRotation', 45)
        title('post')
        legend('topset', 'foreset')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        drawnow
        
    % lineplot of where deposition, normalized to input!
    g = figure();
    subplot(2, 2, 4)
    hold on
        plot(q.Qedg(1:end-1), q.delta_vol.allt' ./ q.qs_inputs, 'o-', 'MarkerSize', 3)
        plot(q.Qedg(1:end-1), q.qs_radflux ./ q.qs_inputs, '^-', 'MarkerSize', 3)
%         plot([1300 1300], [0 1], 'k--')
        plot([min(q.Qedg), max(q.Qedg)], [s.f.lobe_frac, s.f.lobe_frac], 'k--')
        xlabel('water discharge (m^3/s)')
        ylabel('deposited vol. / input vol.')
        legend('delta', 'lobe')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(g, 'Position', [50 100 650 650])
        set(g, 'PaperPositionMode', 'auto')
%         plot(q.Qedg(1:end-1), [q.lobe_vol.preP, q.lobe_vol.postP], '^-', 'MarkerSize', 3) % plot the before prog and after prog breakdown
%         plot(q.Qedg(1:end-1), q.tot_vol.pre ./ sum(q.tot_vol.pre), '+-', 'MarkerSize', 3)
%         plot(q.Qedg(1:end-1), q.tot_vol.post ./ sum(q.tot_vol.post), '+-', 'MarkerSize', 3)
        drawnow; pause(1)
        print('-depsc', '-painters', '-r300',  '../figs/qingshuigou_flux_partitioning.eps');

    % lineplot of where deposition
    g = figure();
    subplot(2, 2, 4)
    hold on
        plot(q.Qedg(1:end-1), [q.delta_vol.alltP], 'o-', 'MarkerSize', 3)
        plot(q.Qedg(1:end-1), [q.lobe_vol.alltP], '^-', 'MarkerSize', 3)
        plot(q.Qedg(1:end-1), q.tot_vol.allt ./ sum(q.tot_vol.allt), '+-', 'MarkerSize', 3)
%         plot([1300 1300], [0 1], 'k--')
        xlabel('water discharge (m^3/s)')
        ylabel('prop. of sediment')
        legend('delta', 'lobe')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(g, 'Position', [50 100 650 650])
        set(g, 'PaperPositionMode', 'auto')
%         plot(q.Qedg(1:end-1), [q.lobe_vol.preP, q.lobe_vol.postP], '^-', 'MarkerSize', 3) % plot the before prog and after prog breakdown
%         plot(q.Qedg(1:end-1), q.tot_vol.pre ./ sum(q.tot_vol.pre), '+-', 'MarkerSize', 3)
%         plot(q.Qedg(1:end-1), q.tot_vol.post ./ sum(q.tot_vol.post), '+-', 'MarkerSize', 3)
%         print('-depsc', '-painters', '-r300',  '../figs/qingshuigou_flux_partitioning.eps');
        
    % lineplot of total volume per bin
    h = figure();
    subplot(2, 2, 4)
        plot(q.Qedg(1:end-1), q.tot_vol.allt./1e9, '+-', 'MarkerSize', 3)
        xlim([0 6e3])
        xlabel('water discharge (m^3/s)')
        ylabel('vol of sed km3')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        set(h, 'Position', [100 100 600 700])
        set(h, 'PaperPositionMode', 'auto')
    
    % calculate the fraction input water at <1300 m3/s (Qwform)
    Qwthresh = 2000;
    Qfidx = q.Qedg <= Qwthresh; % (Qwform)
    lt = sum(q.tot_vol.allt(Qfidx(1:end-1)));
    gt = sum(q.tot_vol.allt(~Qfidx(1:end-1)));
    fr = lt/(lt+gt);
    qwFr = sum(s.Qw<1300) / length(s.Qw);
    disp(['fraction of water discharges less than Qform: ' num2str(qwFr)])
%     disp(['fraction of sediment water discharges less than Qform: ' num2str(fr)])
    
    % do with sediment input
    Qfidx2 = q.Qedg < Qwthresh; % (Qwform)
    lt2 = sum(q.qs_inputs(Qfidx2(1:end-1)));
    gt2 = sum(q.qs_inputs(~Qfidx2(1:end-1)));
    fr2 = lt2/(lt2+gt2);
    qsFr = sum( s.qs(1,s.Qw<1300) ) / sum( s.qs(1,:) );
    disp(['fraction of sediment input for discharges less than Qform: ' num2str(qsFr)])
    
    % display fraction in the lobe
    disp(['fraction of sediment in the lobe at end of run: ' num2str(s.f.lobe_frac)])

   
    % barplot of fraction radflux/influx
    tickidxs = 1:5:length(q.Qcnt);
    f = figure();
    subplot(3, 1, 1)
        bar(1:length(q.Qcnt), [q.qs_inputs', q.qs_radflux'], 1)
%         ylim([0 1])
        set(gca, 'XTick', tickidxs)
        set(gca, 'xTickLabels', q.Qedgelabs(tickidxs), 'xTickLabelRotation', 45)
        title('fluxes')
        legend('input', 'to lobe')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
    subplot(3, 1, 2)
        plot(s.f.tolobe_timeseries ./ s.f.input_timeseries)
        ylim([0 1])
%         set(gca, 'XTick', tickidxs)
%         set(gca, 'xTickLabels', q.Qedgelabs(tickidxs), 'xTickLabelRotation', 45)
        title('fluxes')
%         legend('input', 'to lobe')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
    subplot(3, 1, 3)
        bar(1:length(q.Qcnt), q.qs_radflux' ./ q.qs_inputs', 1)
        ylim([0 1])
        set(gca, 'XTick', tickidxs)
        set(gca, 'xTickLabels', q.Qedgelabs(tickidxs), 'xTickLabelRotation', 45)
        title('fluxes')
%         legend('input', 'to lobe')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        drawnow
    set(f, 'Position', [500 500 450 250])
    set(f, 'PaperPositionMode', 'auto')

end


function [labs] = makelabels(edg)
    labs = cell(1, length(edg)-1);
    for i = 1:length(edg)-2
        labs(i) = {[num2str(edg(i)) '--' num2str(edg(i+1))]};
    end
    labs(end) = {['>' num2str(edg(end-1))]};
end


function [s] = volume_change_per_section(s, d, sig, Bc0, Bf0, Bo0, con)
    
    % need to validate that the depositional width does not change, or
    % calculate incrementally?
    rad_idx = get_idx(d.x, s.rad(1));
    mou_idx = s.mou_idx(end);
    dx = d.dx;
    Be = s.Be(:,end)';
    Bc = s.Bc(:,end)';
    
    lastavul_idx = 1;
    
     % change in eta since last avulsion
%     deta_lastavul = ( s.eta(:, end) + (  s.eta(:, lastavul_idx) - (sig/1000/365*size(s.eta,2))  ) )';
    
    % volume in delta topset region including in channel
    % Vol_inside_rad = ( deta_lastavul(1:rad_idx) .* [dx/2 repmat(dx, 1, rad_idx-2) dx/2] .* Be(1:rad_idx) ) .* (1-con.phi);
    Vol_inside_rad = s.v.Vol_inside_rad;
    
    % volume in the floodplain for redist
    Vol_floodplain = s.v.Vol_inside_floodplain; 
    % ( deta_lastavul(1:rad_idx) .* [dx/2 repmat(dx, 1, rad_idx-2) dx/2] .* (Be(1:rad_idx)-Bc(1:rad_idx)) ) .* (1-con.phi);
    
    detadt = horzcat( zeros(size(s.eta(:,1))), ((s.eta(:, 2:end)+sig*86400) - (s.eta(:, 1:end-1))) );
%     timevar_lobe = sum(detadt(rad_idx:end, :) .* s.Be(rad_idx:end, :), 2);
%     timevar_Vol_lobe_x = timevar_lobe' .* [dx/2 repmat(dx, 1, length(detadt(rad_idx:end, 1))-2) dx/2] * (1-con.phi);
%     timevar_Vol_lobe = sum(timevar_Vol_lobe_x);
%     Vol_lobe = timevar_Vol_lobe;
    Vol_lobe = s.v.Vol_lobe;
    
    % integration approach for the input
    % input_volume = sum(s.qs(1,:)*400 * 86400) * 1e-9; % in km3
    % input_volume_int = sum(trapz(s.qs(1,:)*400) * 86400) * 1e-9;
    input_timeseries = s.qs(1,:) .* s.Bc(1,:) * 86400;
    input_volume = s.v.integrated_input;
    
    % integration approach for the to lobe
    % tolobe_volume = sum(s.qs(rad_idx,:) .* s.Bc(rad_idx,:) * 86400) * 1e-9; % in km3
    % tolobe_volume_int = sum(trapz(s.qs(rad_idx,:) .* s.Bc(rad_idx,:)) * 86400) * 1e-9;
    tolobe_timeseries = s.qs(rad_idx,:) .* s.Bc(rad_idx,:) * 86400;
    tolobe_volume = s.v.Vol_lobe;
    
    % calculate for the whole delta from detadt
    timevar_all_b = (detadt(:, :) .* s.Be(:, :));
    for t = 1:size(timevar_all_b, 2)
        timevar_Vol_all_x(:,t) = timevar_all_b(:,t) .* [0; repmat(dx, length(detadt(:, 1))-1, 1)] .* (1-con.phi);
    end
    s.f.Vdtdx = timevar_Vol_all_x; % Volume per dt per dx
    s.f.Vdtdx_dep = timevar_Vol_all_x; %  make a copy to modify
    s.f.Vdtdx_dep(s.f.Vdtdx_dep < 0) = 0; % Volume per dt per dx only deposition
    
    lobe_cumkm3 = tolobe_volume(end) .* 1e-9;
    top_cumkm3 = Vol_inside_rad(end) .* 1e-9;
    
    lobe_cumkm3perarea = top_cumkm3 / (d.dx * (mou_idx-rad_idx) / 1000);
    top_cumkm3perarea = top_cumkm3 / (d.dx * rad_idx / 1000);

    lobe_cumkm3perlen = top_cumkm3 / (d.dx * (mou_idx-rad_idx) / 1000) / (Bc0+Bo0);
    top_cumkm3perlen = top_cumkm3 / (d.dx * rad_idx / 1000) / (Bc0+Bf0);
    
    tot_cumkm3 = lobe_cumkm3 + top_cumkm3;
    tot_cumkm3perarea = lobe_cumkm3perarea + top_cumkm3perarea;
    tot_cumkm3perlen = lobe_cumkm3perlen + top_cumkm3perlen;
    
    lobe_fracperarea = lobe_cumkm3perarea / tot_cumkm3perarea;
    lobe_fracperlen = lobe_cumkm3perlen / tot_cumkm3perlen;
    lobe_frac = lobe_cumkm3 / tot_cumkm3;
    
    s.f.lobe_cumkm3 = lobe_cumkm3;
    s.f.top_cumkm3 = top_cumkm3;
    
    s.f.lobe_frac = lobe_fracperarea;
    
    s.f.input_timeseries = input_timeseries;
    s.f.tolobe_timeseries = tolobe_timeseries;
   

    
    
end
