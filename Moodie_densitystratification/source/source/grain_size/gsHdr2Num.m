function [gsClassNum] = gsHdr2Num(gsClassHdr)
    % convert header of table of grain size bins to numeric vectors
    class_raw = cellfun(@(x) str2num(x), gsClassHdr, 'Unif', 0); % convert header to numeric, isolates numbers
    class_log = ~cellfun(@(x) isempty(x), class_raw); % logical index for which vars are class bins
    gsClassNum = [class_raw{ class_log }]'; % numeric vector of class edges
end