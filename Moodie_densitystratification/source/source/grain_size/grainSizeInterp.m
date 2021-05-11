function [summNum] = grainSizeInterp(gsDistNum, gsClassNum, queryPts)
% subroutine to makeSummTable     (allows for direct use of this function too)
    
    nullClasses = isnan(gsDistNum); % classes to exclude from interpolation
    zeroClasses = gsDistNum == 0;
    firstFull = find(~or(nullClasses, zeroClasses), 1, 'first');
    if firstFull>1
        gsDistNum(firstFull-1) = 0; % replace one before firstFull with 0
    end
    [summNum] = interp1(cumsum(gsDistNum(~zeroClasses), 'omitnan'), ...
        gsClassNum(~zeroClasses), queryPts); % interpolate for grain sizes at percentiles


end