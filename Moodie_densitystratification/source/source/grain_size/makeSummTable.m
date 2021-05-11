function [summTable, summNum] = makeSummTable(gsDist, queryPts)
    [gsClassNum] = cellfun(@(x) str2double(x), gsDist.Properties.RowNames); % grain size bins numeric
    [gsDistNum] = table2array(gsDist); % grain size dist numeric
    [summNum] = grainSizeInterp(gsDistNum, gsClassNum, queryPts*100); % interpolate for grain sizes at percentiles
    queryStr = join( [ cellstr(repmat('d', length(queryPts), 1)) ...
        cellfun(@(x) num2str(x), num2cell(queryPts * 100), 'Unif', 0)' ], "" ); % unreadable anon functions to make string vector
    [~, lam] = makeDistToData(gsDist);
    summTable = array2table([summNum lam], 'VariableNames', [queryStr; {'lambda'}]); % convert to table
end