function [varMat] = varVsDepth(df, varName, depthDims)
    varCellAll = arrayfun(@(s) eval(['s.' varName]), df, 'Unif', 0);
    varMat = varRepLong(cell2mat(varCellAll), depthDims); % repeat to matrix of same size as concentrations (v = values, m = mat for size)
end