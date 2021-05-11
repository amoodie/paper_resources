function [stationTable] = process_WP_excel_to_struct(river_sheet, key)

    river_filename = ['datasrc/scott_wright_dissertation_data/', river_sheet];

    [status, rawsheets] = xlsfinfo(river_filename);

    datasheets = rawsheets;
    datasheets(ismember( datasheets, {'Summary', 'Chart1'}) ) = [];
    
    stationTable = [];
    for i = 1:length(datasheets)
        disp(['starting: ', num2str(i)])
        
        % main table of stratified
        stationData = readtable(river_filename, 'Sheet', datasheets{i}, 'Range', 'A21:AA171');
        if ismember('Kmdy', stationData.Properties.VariableNames)
            stationData.Kmdy(:) = NaN;
            stationData.Properties.VariableNames{'Kmdy'} = 'Kstdy';
        end
        stationData = table2struct(stationData, 'ToScalar', true);
        
        % other vars table
        otherData = table2array( readtable(river_filename, 'Sheet', datasheets{i}, 'Range', 'B5:B18', 'ReadVariableNames', false) );
        
        stationData.H = otherData(2) * 0.3048; % ft to m
        stationData.ustar = otherData(7);
        stationData.slope = otherData(6);
        stationData.cb = {[stationData.C(1), stationData.C_1(1), stationData.C_2(1)]};
        stationData.D_dim = {[25, 125, 625] .* 1e-6}; %{[88, 177, 354] .* 1e-6}; % {[25, 125, 625] .* 1e-6};
        if ismember(river_sheet, {'MissTar all.xls', 'red all.xls'})
            stationData.d90 = table2array( readtable(river_filename, 'Sheet', datasheets{i}, 'Range', 'J12:J12', 'ReadVariableNames', false) );
        elseif strcmp(river_sheet, 'Niobrara.xls')
                stationData.d90 = otherData(12);
        elseif strcmp(river_sheet, 'RiograndeA2 all.xls')
            stationData.d90 = table2array( readtable(river_filename, 'Sheet', datasheets{i}, 'Range', 'I12:I12', 'ReadVariableNames', false) );
        elseif strcmp(river_sheet, 'midloupE all.xls')
            stationData.d90 = table2array( readtable(river_filename, 'Sheet', datasheets{i}, 'Range', 'F11:F11', 'ReadVariableNames', false) );
        elseif strcmp(river_sheet, 'atchafalaya all.xls')
            stationData.d90 = 0.5;
        else
            error('bad name')
        end
        stationData.ID = datasheets{i};
        stationTable = vertcat(stationTable, stationData);
    end

    river = cell(size(stationTable, 1), 1);
    river(:) = {key};
    [stationTable.river] = river{:};
    
    
end