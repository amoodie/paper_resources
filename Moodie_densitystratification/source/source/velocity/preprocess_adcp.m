function [data] = preprocess_adcp(data)
    % check that the depth range is the same for all data
    try 
        any(~ismember(data.depth(:, 2:end), data.depth(1, 2:end), 'rows'));
    catch
        error('not all depth bins are equal')
    end
    
    % use only columns with more than 50 vel measurements
    nLim = 50;
    [ensemble_idx] = find(data.data(1,1) == data.data(1,:));
    vel_rng = [ensemble_idx(2)+1 ensemble_idx(3)-1];
    [bin_idx] = find(sum(~isnan(data.data(:, vel_rng(1):vel_rng(2)))) > nLim);
    data.data = horzcat(data.data(:, 1:vel_rng(1)-1), ... % metainfo
        data.data(:, bin_idx - 1 + vel_rng(1)), ... % velocity
        data.data(:, vel_rng(2)+1), ... % metainfo
        data.data(:, bin_idx + vel_rng(2) + 1)); % direction
    data.depth = data.depth(1, bin_idx + 1) + data.tddepth; % depth
    data.nbins = length(bin_idx);
    
    % remove adcp velocity outliers -- TRY IT BY A QUALITY INDICATOR FROM ADCP (QA or intensity?)
    velmat = data.data(:, vel_rng(1):vel_rng(1)+data.nbins-1);
    data.IQR = iqr(velmat);
    avg = mean(velmat, 'omitnan');
    low = quantile(velmat, 0.25) - 1.5*data.IQR;
    high = quantile(velmat, 0.75) + 1.5*data.IQR;
    velmat(bsxfun(@lt, velmat, low)) = NaN;
    velmat(bsxfun(@gt, velmat, high)) = NaN;
    data.data(:, vel_rng(1):vel_rng(1)+data.nbins-1) = velmat; % TURN ON TO REMOVE OUTLIERS
end