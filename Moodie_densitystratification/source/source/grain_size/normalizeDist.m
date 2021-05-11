function [gsDistNorm] = normalizeDist(gsDist)
    [newSum] = nansum(gsDist, 1); % new total % of distribution (probably != 100)
    gsDistNorm = gsDist ./ newSum .* 100; % renormalize in interval
end