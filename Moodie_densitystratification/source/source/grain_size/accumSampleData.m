function [sampleStruct] = accumSampleData(sampleStruct, summTable, dataTable, runIdx, gsSummQueryPts)
    sampleStruct.gsSummRawF = summTable(runIdx, :); % summary (with wash load) grain size table all Malvern runs
    sampleStruct.gsDistRawF = dataTable(:, runIdx); % full (with wash load) grain size distribution of all Malvern runs
    [sampleStruct.gsDistMeanF] = averageRuns(sampleStruct.gsDistRawF); % mean and renormalized full grain size distribution
    [sampleStruct.gsDistMeanFNum] = table2array(sampleStruct.gsDistMeanF); % numeric version for compatability issues
    [sampleStruct.gsSummMeanF] = makeSummTable(sampleStruct.gsDistMeanF, gsSummQueryPts); % make summary table of percentile values
end