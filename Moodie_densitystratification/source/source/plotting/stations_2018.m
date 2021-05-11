function stations_2018()

    close all;
    
    load('./dataExport/2018_stations.mat', 'df18')
    df = df18; clear df18
    
    con = load_conset('quartz-water');
    
    global colorSet col1 col2 col3
    colswitch = 'colorful';
    [colorSet, col1, col2, col3] = load_colorSet(colswitch);
    printOut = false;   
    
     %% plot time series data
    stationTime = arrayfun(@(x) x.CollectionDatetime, df);
    ts = readtable('./dataSrc/lijin_discharge_2018.csv');
    % make continuous versions of the discharge and conc vectors
    ts.dischargeCont = ts.discharge;
    nanx = isnan(ts.dischargeCont);
    t = 1:numel(ts.dischargeCont);
    ts.dischargeCont(nanx) = interp1(t(~nanx), ts.dischargeCont(~nanx), t(nanx));
    ts.concentrationCont = ts.concentration;
    nanx = isnan(ts.concentrationCont);
    t = 1:numel(ts.concentrationCont);
    ts.concentrationCont(nanx) = interp1(t(~nanx), ts.concentrationCont(~nanx), t(nanx));
    % mud arrival time
    mudArrivalTime = datetime('07/11/18', 'InputFormat', 'MM/dd/yy');
    
    figure(); hold on
    subplot(2, 1, 1); hold on;
        plot(ts.time, ts.discharge, 'k.', 'MarkerSize', 10)
        plot(ts.time, ts.dischargeCont, 'k-')
        plot(stationTime, interp1(ts.time, ts.dischargeCont, stationTime), ...
            'MarkerFaceColor', col3, 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'LineStyle', 'none')
        plot([mudArrivalTime mudArrivalTime], [1000 4000], 'k--')
        xlim([min(stationTime)-days(3) max(stationTime)+days(3)])
        % ylim([0 100])
        ylabel('discharge (m^3/s)')
        box on;
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')

    subplot(2, 1, 2); hold on;
        plot(ts.time, ts.concentration, 'k.', 'MarkerSize', 10)
        plot(ts.time, ts.concentrationCont, 'k-')
        plot(stationTime, interp1(ts.time, ts.concentrationCont, stationTime), ...
            'MarkerFaceColor', col3, 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'LineStyle', 'none')
        plot([mudArrivalTime mudArrivalTime], [0 30], 'k--')
        xlim([min(stationTime)-days(3) max(stationTime)+days(3)])
        ylabel('concentration (kg/m^3)')
        box on;
        set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')

    %% plot each station by itself
    for i = 1:size(df, 1)
        plot_single_station(df(i))
        
    end
    
   
    
    
end