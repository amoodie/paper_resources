function fit_alpha_function_i()

    source_path = genpath('.');
    data_path = genpath('dataSrc');
    export_path = genpath('dataExport');
    addpath(source_path, data_path, export_path);

    load('toSelectedPlots.mat')
    [con] = load_conset('quartz-water');
    
   
    %% regressions
    figure();
    fmin_func = (@(x) find_alpha(x, reg.Zs, reg.concs, reg.ustarwss, reg.cbis, reg.S, reg.cbSs));
    x0 = [1, 0];
    y = fminsearch(fmin_func, x0);
    
    a = y(1);
    b = y(2);
    
    reg.gamma = a .* (reg.ustarwss) + b;
    
    reg.alphas = reg.gamma .* calculate_alpha_WP04(sum(reg.cbis, 2), reg.S);
    reg.rouses = (1./reg.ustarwss) .* (1./(0.41 .* reg.alphas));
    reg.pred_concs = reg.cbis .* ( ((1-reg.Zs)./reg.Zs).*(0.05./(1-0.05)) ) .^ reg.rouses;
    
    figure()
    plot(reg.pred_concs(:), reg.concs(:), '.')
    
    figure(); 
    subplot(1, 2, 1); hold on
        h1 = plot(reg.concs(:,:,10), reg.Zs);
        set(h1, {'color'}, num2cell(parula(6),2));
        h2 = plot(reg.pred_concs(:,:,10), reg.Zs, '--');
        set(h2, {'color'}, num2cell(parula(6),2));
        h3 = plot(reg.pred_concs0(:,:,10), reg.Zs, ':');
        set(h3, {'color'}, num2cell(parula(6),2));
        set(gca, 'xscale', 'log')
    subplot(1, 2, 2); hold on
        plot(sum(creg.oncs(:,:,10),2), reg.Zs);
        plot(sum(reg.pred_concs(:,:,10),2), reg.Zs, '--');
        plot(sum(reg.pred_concs0(:,:,10),2), reg.Zs, ':');
        colorbar
    
    fig = figure(); hold on;
    plotYRdata(reg.alphas0(:), reg.alphas(:), varRepLong(idx.stationYr(idx.outlierStation), 6), true(282,1), fig, ...
        'colorbyvariable', repmat((1:6)',47,1))
    plot([0 1], [0 1], 'k-')
    colorbar
    colormap(parula(6))
    xlabel('$\alpha_{WP04}$')
    ylabel('$\alpha_{WP04} \cdot \Gamma$')
    
    
    fig = figure(); hold on;
        plotYRdata(reg.ustarwss(:), reg.gamma(:), varRepLong(idx.stationYr(idx.outlierStation), 6), true(282,1), fig, ...
            'colorbyvariable', repmat((1:6)',47,1))
        plot([0 1], [0 1], 'k-')
        colorbar
        colormap(parula(6))
        xlabel('$u_*/w_{si}$')
        ylabel('$\Gamma$')
        
end

function [SSE] = find_alpha(x, Zs, concs, ustarwss, cbis, S, cbSs)

    a = x(1);
    b = x(2);
    
    gamma = a .* ((ustarwss)) + b;
    alphas = gamma .* calculate_alpha_WP04(sum(cbis, 2), S);
    
    rouses = (1./ustarwss) .* (1./(0.41 .* alphas));
    
    pred_concs = cbis .* ( ((1-Zs)./Zs).*(0.05./(1-0.05)) ) .^ rouses;
    
    y = concs(:);
    yhat = pred_concs(:);
    r = y - yhat;
    SSE = ( sum( (yhat-y).^2 ) / length(yhat) ).^(0.5);
    
end

function iter_plot(pred_concs, concs)
    f = gcf();
    hold on;
    cla
    scatter(reshape(pred_concs, [1, 50*47*6]), reshape(concs, [1, 50*47*6]), 55, repelem((1:6)', 50*47, 1), '.')
    plot([1e-16, 1e-0], [1e-16, 1e-0], 'k-')
    set(gca, 'xscale', 'log', 'yscale', 'log')
    xlabel('pred')
    ylabel('meas')
    ylim([1e-5, 1e-2])
    box on
    axis square
    xlim([1e-5, 1e-2])
    drawnow
end
    