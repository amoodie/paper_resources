function [stationStruct] = make_waterSamplesTable(stationStruct)
    [waterSamplesClass] = gsDistClass2num(stationStruct.Samples(1).gsDistMeanNW);
    isWater = ~isnan(vertcat(stationStruct.Samples.metersAboveBed));
    waterSamplesTable = cell(sum(isWater, 1), 11);
    for i = 1:sum(isWater, 1)
        waterSamplesTable{i, 1} = ( stationStruct.Samples(i).NominalDepth ); % depth
        waterSamplesTable{i, 2} = ( stationStruct.Samples(i).NominalZ ) .* stationStruct.FlowDepthAtCollection; % Z
        waterSamplesTable{i, 3} = ( stationStruct.Samples(i).NominalZ ); % Z/H
        waterSamplesTable{i, 4} = stationStruct.Samples(i).IterationNumber; % iteration
        waterSamplesTable{i, 5} = stationStruct.Samples(i).conc; % raw concentrtation
        waterSamplesTable{i, 6} = stationStruct.Samples(i).gsWLdata.gsWLpercent; % washload percentage
        waterSamplesTable{i, 7} = stationStruct.Samples(i).concNW; % concentration with no washload
        waterSamplesTable{i, 8} = waterSamplesClass; % Class
        waterSamplesTable{i, 9} = stationStruct.Samples(i).gsDistMeanNW; % no wash dist
        waterSamplesTable{i, 10} = stationStruct.Samples(i).gsDistMeanNWnorm; % no wash dist norm
        waterSamplesTable{i, 11} = stationStruct.Samples(i).gsSummMeanNWnorm; % no wash dist summaryTable
            waterSamplesTable{i, 12} = table2array( (stationStruct.Samples(i).gsDistMeanNWnorm) )./100 .* stationStruct.Samples(i).concNW; % conc of each class
    end
    [stationStruct.waterSamplesTable] = cell2table(waterSamplesTable, 'VariableNames', ...
        {'sampleDepth', 'sampleZ', 'sampleZnorm', 'sampleIter', 'conc', 'gsWLpercent', 'concNW', 'gsClass', 'gsDistNW', 'gsDistNWnorm', 'gsSummNWnorm', 'concNWbyClass'});
    
end