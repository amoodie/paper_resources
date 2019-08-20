function [idx] =  get_idx(x, loc)
    % get index of some location 'loc' *downstream* from upper boundary
    diffs = abs(x - loc);
    [~, idx] = min(diffs);
    idx = idx - 1;
end