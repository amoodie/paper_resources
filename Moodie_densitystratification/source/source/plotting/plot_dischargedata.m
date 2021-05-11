function plot_dischargedata()

    close all;
    
    load('./dataExport/stationSurveyDataTable.mat')
    
    con = load_conset('quartz-water');
    
    global colorSet col1 col2 col3
    colswitch = 'colorful';
    [colorSet, col1, col2, col3] = load_colorSet(colswitch);
    printOut = false;   
    
    stationTime = arrayfun(@(x) x.CollectionDatetime, df);

    % load data
    t = readtable('lijin_discharge_all.csv');
    
    % interpolations
    t.dischargeCont = t.discharge;
    nanx = isnan(t.dischargeCont);
    tx = 1:numel(t.dischargeCont);
    t.dischargeCont(nanx) = interp1(tx(~nanx), t.dischargeCont(~nanx), tx(nanx));
    
    t15 = t(t.CollectionDatetime < datetime('01/01/2016', 'InputFormat', 'MM/dd/yyyy'), :);
    t16 = t(and(t.CollectionDatetime > datetime('01/01/2016', 'InputFormat', 'MM/dd/yyyy'), t.CollectionDatetime < datetime('01/01/2017', 'InputFormat', 'MM/dd/yyyy')), :);
    t18 = t(t.CollectionDatetime > datetime('01/01/2018', 'InputFormat', 'MM/dd/yyyy'), :);
    
    figure(); hold on
    subplot(1, 3, 1); hold on;
        plot(t15.CollectionDatetime, t15.discharge, 'k.', 'MarkerSize', 10)
        plot(t15.CollectionDatetime, t15.dischargeCont, 'k-')
        plot(stationTime, interp1(t.CollectionDatetime, t.dischargeCont, stationTime), ...
            'MarkerFaceColor', col1, 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'LineStyle', 'none')
        ylim([0 4000])
        ylabel('discharge (m^3/s)')
        xlim(datetime({'07/02/2015', '07/15/2015'}, 'InputFormat', 'MM/dd/yyyy'))
        box on;
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
    subplot(1, 3, 2); hold on;
        plot(t16.CollectionDatetime, t16.discharge, 'k.', 'MarkerSize', 10)
        plot(t16.CollectionDatetime, t16.dischargeCont, 'k-')
        plot(stationTime, interp1(t.CollectionDatetime, t.dischargeCont, stationTime), ...
            'MarkerFaceColor', col2, 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'LineStyle', 'none')
        ylim([0 4000])
%         ylabel('discharge (m^3/s)')
        xlim(datetime({'07/02/2016', '07/15/2016'}, 'InputFormat', 'MM/dd/yyyy'))
        box on;
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
    subplot(1, 3, 3); hold on;
        plot(t18.CollectionDatetime, t18.discharge, 'k.', 'MarkerSize', 10)
        plot(t18.CollectionDatetime, t18.dischargeCont, 'k-')
        plot(stationTime, interp1(t.CollectionDatetime, t.dischargeCont, stationTime), ...
            'MarkerFaceColor', col3, 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'LineStyle', 'none')
        mudArrivalTime = datetime('07/11/2018', 'InputFormat', 'MM/dd/yyyy');
        plot([mudArrivalTime mudArrivalTime], [0 4000], 'k--')
        xlim(datetime({'07/02/2018', '07/15/2018'}, 'InputFormat', 'MM/dd/yyyy'))
        ylim([0 4000])
%         ylabel('discharge (m^3/s)')
        box on;
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')


end