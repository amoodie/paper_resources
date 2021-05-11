function [stationStruct] = make_simpleDatasets(stationStruct, simplifyNclasses)

    %% temporary store old dists for comparison
    origBedDist = stationStruct.bedData.gsDistBed;
    origNearBedDist = stationStruct.nearBedData.gsDistNearBedNWnorm;
    

    %% water samples

    % make concentration profiles with subsets of the grain size profiles
    %   make a mroe complex rebinning process that can handle arbitrary classes.
    %   determine the bin spacing by: logspace(log10(3),log10(10),10) (ten logspaced bins 3 to 10).
    %   find the bins to index on the fly and grab those, choose the *midpoint* instead of the initial binning and put this in the table row names
    binFlag = 'manual'; % 'linear', 'log', 'manual'
    if strcmp(binFlag, 'linear')
        gsBins = stationStruct.waterSamplesTable.gsClass{1};
        lims = round(linspace(gsBins(2), 350, simplifyNclasses), 3);
        numElems = zeros(1, length(gsBins));
        for i = length(lims):-1:1
            numElems(gsBins <= lims(i)) = i;
        end
    elseif strcmp(binFlag, 'log')
        lims = round(logspace(1.2, 2.69897, simplifyNclasses), 3);
        gsBins = stationStruct.waterSamplesTable.gsClass{1};
        numElems = zeros(1, length(gsBins));
        for i = length(lims):-1:1
            numElems(gsBins <= lims(i)) = i;
        end
    elseif strcmp(binFlag, 'manual')
        lims = [25, 43, 73, 122, 208, 350];
        gsBins = stationStruct.waterSamplesTable.gsClass{1};
        numElems = zeros(1, length(gsBins));
        for i = length(lims):-1:1
            numElems(gsBins <= lims(i)) = i;
        end
    end
    
    % preallocate
    gsClass = NaN(size(stationStruct.waterSamplesTable,1), simplifyNclasses);
    gsChar = NaN(size(stationStruct.waterSamplesTable,1), simplifyNclasses);
    gsDistNW = NaN(size(stationStruct.waterSamplesTable,1), simplifyNclasses);
    gsDistNWnorm = NaN(size(stationStruct.waterSamplesTable,1), simplifyNclasses);
    concNWbyClass = NaN(size(stationStruct.waterSamplesTable,1), simplifyNclasses);

        for s = 1:size(stationStruct.waterSamplesTable, 1)
        % make the reduced bin versions of each element in the water table
            gsDistNWcum(s, :) = interp1(gsDistClass2num(stationStruct.waterSamplesTable.gsDistNW{s}), ...
                                     cumsum(stationStruct.waterSamplesTable.gsDistNW{s}.avgClassNorm, 'omitnan'), lims);
            gsDistNWnormCum(s, :) = interp1(gsDistClass2num(stationStruct.waterSamplesTable.gsDistNWnorm{s}), ...
                                     cumsum(stationStruct.waterSamplesTable.gsDistNWnorm{s}.avgClassNorm, 'omitnan'), lims);
            concNWbyClassCum(s, :) = interp1(gsDistClass2num(stationStruct.waterSamplesTable.gsDistNWnorm{s}), ...
                                     cumsum(stationStruct.waterSamplesTable.concNWbyClass{s}, 'omitnan'), lims);
            gsClass(s, :) = lims;
            
        end
        gsDistNW = [gsDistNWcum(:,1), gsDistNWcum(:,2:end) - gsDistNWcum(:,1:end-1)];
        gsDistNWnorm = [gsDistNWnormCum(:,1), gsDistNWnormCum(:,2:end) - gsDistNWnormCum(:,1:end-1)];
        concNWbyClass = [concNWbyClassCum(:,1), concNWbyClassCum(:,2:end) - concNWbyClassCum(:,1:end-1)];
        gsChar = gsClass;

    gsDistNearBedNWcum = interp1(gsDistClass2num(stationStruct.nearBedData.gsDistNearBedNW), ...
                                     cumsum(stationStruct.nearBedData.gsDistNearBedNW.avgClassNorm, 'omitnan'), lims);
    gsDistNearBedNWnormCum = interp1(gsDistClass2num(stationStruct.nearBedData.gsDistNearBedNWnorm), ...
                                     cumsum(stationStruct.nearBedData.gsDistNearBedNWnorm.avgClassNorm, 'omitnan'), lims);
    gsDistNearBedNW = [gsDistNearBedNWcum(1), gsDistNearBedNWcum(2:end) - gsDistNearBedNWcum(1:end-1)];
    gsDistNearBedNWnorm = [gsDistNearBedNWnormCum(1), gsDistNearBedNWnormCum(2:end) - gsDistNearBedNWnormCum(1:end-1)];
    
    % loop through each sample, replacing with relevent reduced bins
    for s = 1:size(stationStruct.waterSamplesTable, 1)
        
        % gsClass
        stationStruct.waterSamplesTable.gsClass(s) = { gsClass(s, :)' };
        
        % gsChar
        stationStruct.waterSamplesTable.gsChar(s) = { gsChar(s, :)' };
        
        % gsDistNW
        stationStruct.waterSamplesTable.gsDistNW{s} = array2table(gsDistNW(s,:)', 'RowNames', strtrim(cellstr(num2str(gsClass(s,:)'))), ...
            'VariableNames', stationStruct.waterSamplesTable.gsDistNW{s}.Properties.VariableNames(1)); % convert to table with var names
        
        %gsDistNWnorm
        stationStruct.waterSamplesTable.gsDistNWnorm{s} = array2table(gsDistNWnorm(s,:)', 'RowNames', strtrim(cellstr(num2str(gsClass(s,:)'))), ...
            'VariableNames', stationStruct.waterSamplesTable.gsDistNWnorm{s}.Properties.VariableNames(1)); % convert to table with var names
        
        % concNWbyClass
        stationStruct.waterSamplesTable.concNWbyClass{s} =  concNWbyClass(s,:)' ;
        
    end
    
    % process to the nearBedData struct
    stationStruct.nearBedData.gsDistNearBedNW = array2table(gsDistNearBedNW', 'RowNames', strtrim(cellstr(num2str(gsClass(s,:)'))), ...
        'VariableNames', stationStruct.nearBedData.gsDistNearBedNW.Properties.VariableNames(1)); % convert to table with var names
    stationStruct.nearBedData.gsDistNearBedNWnorm = array2table(gsDistNearBedNWnorm', 'RowNames', strtrim(cellstr(num2str(gsClass(s,:)'))), ...
        'VariableNames', stationStruct.nearBedData.gsDistNearBedNWnorm.Properties.VariableNames(1)); % convert to table with var names
    stationStruct.nearBedData.gsCharNearBedNWnorm = array2table(gsDistNearBedNWnorm', 'RowNames', strtrim(cellstr(num2str(gsChar(s,:)'))), ...
        'VariableNames', stationStruct.nearBedData.gsDistNearBedNWnorm.Properties.VariableNames(1)); % convert to table with var names
    stationStruct.nearBedData.gsSummNearBedNWnorm = stationStruct.nearBedData.gsSummNearBedNWnorm;
    
    
    %% bed samples
    
    % kill zeros at front of distribution
    firstnonzero = find(stationStruct.bedData.gsDistBed.avgClassNorm > 0, 1, 'first');
    stationStruct.bedData.gsDistBed = stationStruct.bedData.gsDistBed(firstnonzero:end, :);
    
    % make concentration profiles with subsets of the grain size profiles
    % binFlag = 'manual'; % 'linear', 'log', 'manual'
    if strcmp(binFlag, 'linear')
        gsBins = gsDistClass2num(stationStruct.bedData.gsDistBed);
        lims = round(linspace(gsBins(2), 350, simplifyNclasses), 3);
        numElems = zeros(1, length(gsBins));
        for i = length(lims):-1:1
            numElems(gsBins <= lims(i)) = i;
        end
    elseif strcmp(binFlag, 'log')
        lims = round(logspace(1.2, 2.69897, simplifyNclasses), 3);
        gsBins = gsDistClass2num(stationStruct.bedData.gsDistBed);
        lims(1) = max([lims(1), gsBins(1)]);
        numElems = zeros(1, length(gsBins));
        for i = length(lims):-1:1
            numElems(gsBins <= lims(i)) = i;
        end
    elseif strcmp(binFlag, 'manual')
        % must be exactly the same as near bed, in order for entrainment to work
        lims = gsClass(1,:);
        gsBins = gsDistClass2num(stationStruct.bedData.gsDistBed);
        numElems = zeros(1, length(gsBins));
        for i = length(lims):-1:1
            numElems(gsBins <= lims(i)) = i;
        end
    end
    
    bedDistClasses = gsDistClass2num(stationStruct.bedData.gsDistBed);
    bedDistCumDist = cumsum(stationStruct.bedData.gsDistBed.avgClassNorm);
    lastfrac = 0;
    for c = 1:simplifyNclasses
            
        % determine the value to fill for fraction
        if lims(c) < bedDistClasses(1)
            newcum = 0;
        elseif lims(c) > bedDistClasses(end)
            newcum = 100.0;
        else
            newcum = interp1(bedDistClasses, bedDistCumDist, lims(c));
        end
        fractionc = newcum - lastfrac;
        lastfrac = newcum;
        
        [gsClassBed(c, 1)] = lims(c);
        [gsDistBed(c, 1)] = fractionc; 
        
    end
    
    % process to the bedData struct
    stationStruct.bedData.gsDistBed = array2table(gsDistBed, 'RowNames', strtrim(cellstr(num2str(gsClassBed))), ...
        'VariableNames', stationStruct.bedData.gsDistBed.Properties.VariableNames(1)); % convert to table with var names
    stationStruct.bedData.gsCharBed = array2table(gsDistBed, 'RowNames', strtrim(cellstr(num2str(gsClassBed))), ...
        'VariableNames', stationStruct.bedData.gsDistBed.Properties.VariableNames(1)); % convert to table with var names
    
    if false
        figure(); hold on;
        plot(gsDistClass2num(origBedDist), cumsum(origBedDist.avgClassNorm, 'omitnan'), 'ok-')
        plot(gsDistClass2num(stationStruct.bedData.gsDistBed), cumsum(stationStruct.bedData.gsDistBed.avgClassNorm, 'omitnan'), 'or-')

        plot(gsDistClass2num(origNearBedDist), cumsum(origNearBedDist.avgClassNorm, 'omitnan'), 'ok--')
        plot(gsDistClass2num(stationStruct.nearBedData.gsDistNearBedNWnorm), cumsum(stationStruct.nearBedData.gsDistNearBedNWnorm.avgClassNorm, 'omitnan'), 'or--')
    end
        
end