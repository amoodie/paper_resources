function [SampleIDlist, summTab, dataTab, gsClass] = gsRawDataLoad()
    filename = 'grain-size_distributions.csv'; % distribution data of all samples
    tb = readtable(filename, 'delimiter', ',', 'ReadVariableNames', 1); % read in table
    fid = fopen(filename); % open file id
    hdr_raw = textscan(fid, '%s', (size(tb, 2)), 'delimiter', ','); % read header line
    hdr = hdr_raw{:}; % process header line
    class_raw = cellfun(@(x) str2num(x), hdr, 'Unif', 0); % convert header to numeric, isolates numbers
    class_log = ~cellfun(@(x) isempty(x), class_raw); % logical index for which vars are class bins
    gsClass = [class_raw{ class_log }]'; % numeric vector of class edges
    dX_log = cellfun(@(x) x(1), hdr) == 'd'; % logical index for summary table values
    
    SampleIDlist = tb.SampleName; % list of sample names
    summTab = tb(:, dX_log); % list of summary table values
    dataTab = tb(:, class_log); % all distribution data
    dataTab = table2array(dataTab)'; % transpose the data table to have grain size in rows
    dataTab = array2table(dataTab, 'RowNames', strtrim(cellstr(num2str(gsClass))), ...
        'VariableNames', matlab.lang.makeUniqueStrings( matlab.lang.makeValidName(SampleIDlist)) ); % convert back to table with row and var names
end