function [gsDist] = averageRuns(gsDistRaw)
    if ~isempty(gsDistRaw)
        [avgClass] = nanmean(table2array(gsDistRaw), 2); % average % in each gsClass
        [avgClassNorm] = normalizeDist(avgClass); % renormalize the gsDist to 100%
        [gsDist] = array2table(avgClassNorm, 'RowNames', gsDistRaw.Properties.RowNames); % convert to table with var names
    else
        [gsDist] = gsDistRaw;
    end
end