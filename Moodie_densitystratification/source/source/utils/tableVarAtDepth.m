function [depthVar, depthVarDims] = tableVarAtDepth(df, varName, depth)
    %tableVarAtDepth get a variable for all samples at a given depth
    %
    % [depthVar, depthVarDims] = tableVarAtDepth(df, varName, depth)

    depthIdx = arrayfun(@(s) s.waterSamplesTable.sampleDepth == depth, df, 'Unif', 0);
    depthVarAll = arrayfun(@(s) eval(['s.waterSamplesTable.' varName]), df, 'Unif', 0);
    depthVarCell = cellfun(@(s, idx) s(idx), depthVarAll, depthIdx, 'Unif', 0); % extract the near bed concentrations only (s = struct, idx = index)
    depthVar = vertcat(depthVarCell{:}); % reformat into long list for plotting
    depthVarDims = cellfun(@length, depthVarCell); % number of near bed concs from each station

end