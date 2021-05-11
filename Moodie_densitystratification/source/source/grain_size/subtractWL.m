function [sampleStruct] = subtractWL(sampleStruct, gsDistBed, gsPercWLcutoff, gsGSWLcutoff) % grain size distribution without washload
    
    % check that only one of the options is given!
    if all(isnan([gsPercWLcutoff, gsGSWLcutoff]))
        error('must define either gsPercWLcutoff or gsGSWLcutoff, both are undef now')
    end
    if all(~isnan([gsPercWLcutoff, gsGSWLcutoff]))
        error('must define either gsPercWLcutoff or gsGSWLcutoff, both are given now')
    end
    
    % process the option given to a grain size
    if isnan(gsGSWLcutoff)
        [~, gsWLcutoff] = makeSummTable(gsDistBed, gsPercWLcutoff); % get washload cutoff grainsize
        bedSampleClass = gsDistClass2num(gsDistBed);
    else
        gsWLcutoff = gsGSWLcutoff;
    end
    
    waterSampleClass = gsDistClass2num(sampleStruct.gsDistMeanF);
    waterSampleDistNum = table2array(sampleStruct.gsDistMeanF);
    [gsWLpercent] = interp1(waterSampleClass, cumsum(table2array(sampleStruct.gsDistMeanF), 'omitnan'), gsWLcutoff); % percent of dist below washload cutoff
    [gsWLupperIdx] = find(waterSampleClass > gsWLcutoff, 1, 'first'); % index of first bin larger than washload cutoff
    newClass = [gsWLcutoff; waterSampleClass(gsWLupperIdx:end)];
    newDist = [NaN; waterSampleDistNum(gsWLupperIdx:end)];
    newDistToNorm = [NaN; waterSampleDistNum(gsWLupperIdx:end)];
    [sampleStruct.gsDistMeanNW] = array2table(newDist, 'RowNames', strtrim(cellstr(num2str(newClass))), ...
        'VariableNames', sampleStruct.gsDistMeanF.Properties.VariableNames(1)); % convert to table with var names
    [sampleStruct.gsDistMeanNWnorm] = array2table(normalizeDist(newDistToNorm), 'RowNames', strtrim(cellstr(num2str(newClass))), ...
        'VariableNames', sampleStruct.gsDistMeanF.Properties.VariableNames(1)); % convert to table with var names
    [sampleStruct.gsWLdata] = array2table([gsWLcutoff, gsWLpercent], 'VariableNames', {'gsWLcutoff', 'gsWLpercent'});
    concVol = sampleStruct.conc / 2.65 / 1000; % volumetric concentration, needed to use % volume distribution
    concVolNW = concVol * ((100 - gsWLpercent)/100);
    sampleStruct.concNW = concVolNW * 2.65 * 1000;
end