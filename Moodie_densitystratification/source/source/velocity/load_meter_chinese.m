function [meter] = load_meter_chinese(stnInfoMod)
    
    filename = strcat('./datasrc/velocimeter_csvs/', stnInfoMod, '_format.csv');
    
    fid = fopen(filename);
    tline = fgetl(fid);
    fclose(fid);
    ncols = length(find( tline == ',' )) + 1;
    
    meter.data = csvread(filename, 9, 1, [9, 1, 42, ncols-1]);
    meter.internalcalib = 1 ./ (meter.data(2:end, :) - meter.data(1:end-1, :));
    meter.rotpsec = 20 ./ (meter.data(2:end, :) - meter.data(1:end-1, :));
    
    meter.flowDepth = csvread(filename, 2, 1, [2 1 2 1]);
    meter.totalTime = csvread(filename, 45, 1, [45 1 45 ncols-1]);
    meter.velocityFromField = csvread(filename, 46, 1, [46 1 46 ncols-1]);
    
    meter.measDepth = csvread(filename, 7, 1, [7 1 7 ncols-1]);
    meter.measZ = meter.flowDepth - meter.measDepth;
    
    % calibration for high flow from chinese colleagues is used
    calib_from_china = (@(x) 0.2540 .* x + 0.0035);
    % std= +- 0.97%
    
    % this is some other relationship I worked up in the interrum, based on
    % the internal log data
    calib_linear_fit = (@(x) 4.999 .* x + 0.02947);
    calib_linear_fit_ll = (@(x) 4.869 .* x + -0.03983);
    calib_linear_fit_ul = (@(x) 5.129 .* x + 0.09878);
    % 95% CI is: (4.869, 5.129) and (-0.03983, 0.09878)

    meter.veldata = calib_from_china(meter.rotpsec);
    
    % if using the calibration from colleagues
    meter.velmean = nanmean(meter.veldata);
    meter.velmedian = median(meter.veldata, 'omitnan');
    meter.velerr = std(meter.veldata, 'omitnan');
    meter.vellerr_key = 'std';

    if strcmp(stnInfoMod, 'KC1_070918') && false
        m = calib_linear_fit(meter.internalcalib);
        l = calib_linear_fit_ll(meter.internalcalib);
        u = calib_linear_fit_ul(meter.internalcalib);
        lm = calib_linear_fit_ll(nanmean(meter.internalcalib));
        um = calib_linear_fit_ul(nanmean(meter.internalcalib));
        
        figure(); hold on;
        plot(meter.internalcalib, meter.velocityFromField, 'ko')
        plot(meter.internalcalib, u, 'r-')
        plot(meter.internalcalib, l, 'b-')
        plot(meter.internalcalib, m, 'k-', 'LineWidth', 2)
        plot(nanmean(meter.internalcalib), calib_linear_fit_ll(nanmean(meter.internalcalib)), 'b.', 'MarkerSize', 12)
        plot(nanmean(meter.internalcalib), calib_linear_fit_ul(nanmean(meter.internalcalib)), 'r.', 'MarkerSize', 12)
        plot(nanmean(meter.internalcalib), calib_linear_fit(nanmean(meter.internalcalib)), 'k.', 'MarkerSize', 30)
        ylabel('velocity (m/s)')
        xlabel('n per second (20 rotations/s)')
        title('station KC1\_070918')
    end
    
end