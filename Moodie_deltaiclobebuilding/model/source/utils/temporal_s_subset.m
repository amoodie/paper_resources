function [s] = temporal_s_subset(s, idx)
%% temporal subset of structure s, idx is logical of time length of s

    names = fieldnames(s);
    for i = 1:numel(names)
        s.(names{i}) = subset(s.(names{i}), idx);
    end

    end

function [ss] = subset(x, idx)
    if isvector(x)
        ss = x(idx);
    else
        ss = x(:, idx);
    end
end