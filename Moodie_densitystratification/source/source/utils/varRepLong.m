function [repLong] = varRepLong(varMat, dims)
    %varRepLong repeats a variable to match the number of dims needed for multi-sample at depth analysis
    %
    % [repLong] = varRepLong(varMat, dims)
    
    if isscalar(dims)
        dims = repmat(dims, size(varMat, 1), 1);
    end

    repLongCell = arrayfun(@(v, r) repmat(v, r, 1), varMat, dims, 'Unif', 0);
    repLong = vertcat( repLongCell{:} );
end