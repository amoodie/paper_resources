function [gsClassNum] = gsDistClass2num(gsClassTable)
    %gsDistClass2num convert rownames to numeric

    [gsClassNum] = cellfun(@(x) str2double(x), gsClassTable.Properties.RowNames); % numeric of gsClass values

end

