function [meter] = load_meter_swoffer(stnInfoMod)
    
    filename = strcat('./datasrc/velocimeter_csvs/', stnInfoMod, '_format.csv');
    
    meter.data = csvread(filename, 7, 1);
    meter.veldata = meter.data;
    meter.flowDepth = csvread(filename, 2, 1, [2 1 2 1]);
    
    N = size(meter.data, 2);
    meter.measZ = csvread(filename, 5, 1, [5 1 5 N]);
    
    meter.velmean = mean(meter.veldata, 'omitnan');
    meter.velmedian = median(meter.veldata, 'omitnan');
    meter.velerr = std(meter.veldata, 'omitnan');
    meter.vellerr_key = 'std';
    
end