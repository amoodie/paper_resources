function [adcp] = process_adcp(adcp)

    adcp.mean = mean(adcp.data(:, 6:6+adcp.nbins-1), 'omitnan');
    adcp.median = median(adcp.data(:, 6:6+adcp.nbins-1), 'omitnan');
    adcp.std = std(adcp.data(:, 6:6+adcp.nbins-1), 'omitnan');

end