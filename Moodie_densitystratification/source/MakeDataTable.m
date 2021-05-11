function [df] = MakeDataTable(varargin)
% MakeDataTable    master file
% This is the master file for processing all data from 2015 and 2016 and 2018 surveys
    % into a single large data table. The file first creates the template 
    % for the table, then fills it in with the station data. 
    % 
    % Options available are only to either process or load sections of the 
    % dataset. Default behavior is to process each portion of the dataset 
    % into the table on execution of this file.
    % 
    % 
    source_path = genpath('source');
    data_path = genpath('datasrc');
    addpath(source_path, data_path);

    %% Parse the input variables (if present)
    % check inputs and clean
    if nargin > 0
        typechk = all(ismember(horzcat(varargin{2:2:end}), [0 1])); % are all binary 0/1
        if ~typechk
            error('input values not of type logical or convertable to logical');
        end
        varargin(2:2:end) = num2cell(  logical( horzcat(varargin{2:2:end}) )  ); % convert varargins to logical
    end
    
    % define and execute parse
    defaultOption = true;
    p = inputParser;
    addParameter(p, 'stnInfo', defaultOption);
    addParameter(p, 'velocity', defaultOption);
    addParameter(p, 'grainSize', defaultOption);
    addParameter(p, 'concentration', defaultOption);
    
    parse(p, varargin{:})
    
    
    % define additional parameters for how to process the data
    simplifyFlag = true; % simplify the grain size distributions into fewer classes
    simplifyNclasses = 6;
    
    globalAverageBedDistribution = false; % replace all bed material samples with a global average for smoothing
    
    gsSummQueryPts = [0.1 0.5 0.9 0.84 0.05 0.95]; % summary table query points
    gsPercWLcutoff = NaN;% 0.02; % 0.05; % washload cutoff percentile
    gsGSWLcutoff = 15; % washload cutoff grain size
    % note that an error is thrown if both gsPercWLcutoff and gsGSWLcutoff
    % are given, one must be NaN!!

    warning('off','MATLAB:table:ModifiedAndSavedVarnames')
    warning('off','MATLAB:textscan:UnableToGuessFormat')

    
    %% stnInfo
    % process/load
    % binary for processing (1) or loading (0)
    if p.Results.stnInfo
        tb = readtable('concentration_measurements.csv', 'delimiter', ',', 'DateLocale', 'en_US');
        [stnIDs, stnIDsIdx] = unique(tb.StationID);
        stnInfoExtractVars = {'StationID', 'NominalLocation', 'Transect', 'WrtChannel', ...
            'StationLocation', 'CollectionDate', 'CollectionTime', 'FlowDepthAtCollection'};
        for i = 1:length(stnIDs)
            stnInfoRaw(i, :) = tb(stnIDsIdx(i), stnInfoExtractVars);
        end
        df = table2struct(stnInfoRaw);
        CollectionDatetime = arrayfun(@(x) datetime(datenum([x.CollectionDate ' ' x.CollectionTime]), 'ConvertFrom', 'datenum'), df, 'Unif', 0);
        [df.CollectionDatetime] = CollectionDatetime{:};
    else 
        % loading function here?
    end
    
    %% concentration
    if p.Results.concentration
        tb = readtable('concentration_measurements.csv', 'delimiter', ',');
        for i = 1:size(df, 1)
            df(i).uncertaintyTable = [];
            matchIdx = strcmp(tb.StationID, tb{stnIDsIdx(i), 'StationID'});
            sampleInfoExtractVars = {'FullSampleName', 'SampleID', 'SampleLocation', 'WrtChannel', ...
                'NominalDepth', 'NominalZ', 'IterationNumber', 'FlowDepthAtCollection', ...
                'PredictedCollectionDepth', 'ActualCollectionDepth', 'SmoothedCollectionDepth', 'metersAboveBed', ...
                'TotalVolume_mL_', 'TotalMass_g_', 'Concentration_g_l_'};
            concDataExtractVars = {'GraduatedCylinder', 'StartTimeOfFilter', 'TotalVolume_mL_', ...
                'EmptyFilterWeight_g_', 'FilterUnitNumber', 'EndTimeOfFilter', 'FullFilterWeight_dry__g_', ...
                'dFilter_g_', 'WeightOfEmptyBeaker_g_', 'WeightOfFullBeaker_dry__g_', 'MassOfSedimentInBeaker_g_', ...
                'TotalMass_g_', 'PercentOnFilter___', 'Concentration_g_l_', 'concErrorDown', 'concErrorUp'};
            sampleInfoRaw = tb(matchIdx, sampleInfoExtractVars);
            concDataRaw = tb(matchIdx, concDataExtractVars);
            % rename a few vars:
                sampleInfoRaw.Properties.VariableNames(strcmp(sampleInfoRaw.Properties.VariableNames, 'TotalVolume_mL')) = {'TotalVolume_mL'};
                sampleInfoRaw.Properties.VariableNames(strcmp(sampleInfoRaw.Properties.VariableNames, 'TotalMass_g_')) = {'TotalMass_g'};
                sampleInfoRaw.Properties.VariableNames(strcmp(sampleInfoRaw.Properties.VariableNames, 'Concentration_g_l_')) = {'RawConcentration_gl'};
                if any(or( ~isnan(concDataRaw.concErrorUp), ~isnan(concDataRaw.concErrorDown) ))
                    sampleInfoRaw.conc = nansum([sampleInfoRaw.RawConcentration_gl, concDataRaw.concErrorUp, -1*concDataRaw.concErrorDown], 2);
                    df(i).uncertaintyTable.sampleConcMod = or( ~isnan(concDataRaw.concErrorUp), ~isnan(concDataRaw.concErrorDown) );
                else
                    sampleInfoRaw.conc = sampleInfoRaw.RawConcentration_gl;
                end
            sampleInfoRaw.concData = (concDataRaw);
            df(i).Samples = table2struct(sampleInfoRaw);
        end
    else
        % loading function here?
    end
    
    %% velocity
    discharge_table = readtable('lijin_discharge_all.csv');
    if p.Results.velocity
        for i = 1:size(df, 1)
            stni = strtrim(tb{stnIDsIdx(i), 'StationID'});
            splitStn = strsplit(stni{:}, '/');
            yrstrraw = splitStn{end};
            yrstr = yrstrraw(~isletter(yrstrraw));
            yr = str2double(yrstr);
            if yr == 16
                [df(i).Velocity] = MakeVelocity16(stni{:}, df(i));
                df(i).Velocity.flowDepth = df(i).Samples(1).FlowDepthAtCollection;
                df(i).Velocity.discharge = load_discharge(df(i).CollectionDatetime, discharge_table);
            elseif yr == 18
                [df(i).Velocity] = MakeVelocity18(stni{:}, df(i));
                df(i).Velocity.flowDepth = df(i).Samples(1).FlowDepthAtCollection;
                df(i).Velocity.discharge = load_discharge(df(i).CollectionDatetime, discharge_table);
            elseif yr == 15
                df(i).Velocity.discharge = load_discharge(df(i).CollectionDatetime, discharge_table);
                df(i).Velocity.meter = [];
                df(i).Velocity.adcp = [];
                df(i).Velocity.data.measZ = [];
                df(i).Velocity.flowDepth = df(i).Samples(1).SmoothedCollectionDepth;
            end
            % do following to all years
            df(i).Velocity.slope = 6.4e-5;
            df(i).Velocity.Cf = 0.001; % Ma et al., 2017
            df(i).Velocity.taubDSP = df(i).Velocity.flowDepth * df(i).Velocity.slope * 9.81 * 1000;
            df(i).Velocity.ustarDSP = sqrt( df(i).Velocity.flowDepth * df(i).Velocity.slope * 9.81 );
            df(i).Velocity.Ubar_crosssection = df(i).Velocity.discharge / (400 * df(i).Velocity.flowDepth);
            load('../data/ustar_calibration.mat', 'ustar_calib_func')
            df(i).Velocity.ustarCalib = ustar_calib_func(df(i).Velocity.discharge, df(i).Velocity.flowDepth);
        end
        % calculate a transect average ustar to smooth things
        [df] = transectAvgUstar(df);
    end
    
    %% grainSize
    meanBedDist = load('meanBedDist.mat');
    meanBedDist = meanBedDist.meanBedDist;
    meanBedDistSumm = makeSummTable( meanBedDist, gsSummQueryPts );
    if p.Results.grainSize
        [sampleIDlist, summTable, dataTable, gsOrigClass] = gsRawDataLoad();
        [sampleIDlistStripped] = cellfun(@(x) x{1}, cellfun(@(x) strsplit(x, '[a-z]', 'DelimiterType', 'RegularExpression'), ...
            sampleIDlist, 'Unif', 0), 'Unif', 0); % split string "iter" off end for no-gs dists
        for i = 1:size(df, 1)
            % preallocate the fields of df as needed
            samplesBed = find(isnan( vertcat(df(i).Samples.metersAboveBed) ));
            samplesWater = find(~isnan( vertcat(df(i).Samples.metersAboveBed) ));
            df(i).Samples(end).gsSummRawF = [];
            df(i).Samples(end).gsDistRawF = [];
            df(i).Samples(end).gsDistMeanF = [];
            df(i).Samples(end).gsDistMeanFNum = [];
            df(i).Samples(end).gsSummMeanF = [];
            df(i).Samples(end).gsWLdata = [];
            df(i).Samples(end).gsDistMeanNW = [];
            df(i).Samples(end).gsDistMeanNWnorm = [];
            df(i).Samples(end).concNW = [];
            df(i).Samples(end).gsSummMeanNWnorm = [];
            df(i).waterSamplesTable = [];
            df(i).nearBedData = [];
            df(i).bedData = [];
            
            % accum bed samples first
            if ~isempty(samplesBed)
                for j = 1:length(samplesBed)
                    sampleIdx = samplesBed(j);
                    sampleID = df(i).Samples(sampleIdx).SampleID;
                    runIdx = strcmp( sampleIDlist, sampleID );
                    if any(runIdx) % if there is grain size data for this sample specifically, skip to processing
                    elseif isnan(df(i).Samples(sampleIdx).metersAboveBed) % there is no grain size data for this sample and this is a bed sample
                        sampleIDStripped = sampleID(1:end-1);
                        runIdx = strcmp( sampleIDlistStripped, sampleIDStripped );
                    end
                    [df(i).Samples(samplesBed)] = accumSampleData(df(i).Samples(sampleIdx), summTable, dataTable, runIdx, gsSummQueryPts);
                    df(i).uncertaintyTable.notActualBed = 0;
                end
                try
                    if globalAverageBedDistribution
                        [df(i).gsDistBed] = meanBedDist;
                        [df(i).gsSummBed] = meanBedDistSumm;
                    else
                        df(i).bedData.gsDistBedF = array2table([df(i).Samples(samplesBed).gsDistMeanFNum], ...
                            'RowNames', df(i).Samples(sampleIdx).gsDistMeanF.Properties.RowNames, ...
                            'VariableNames', matlab.lang.makeValidName(matlab.lang.makeUniqueStrings( {df(i).Samples(samplesBed).SampleID} )) );
                        [df(i).bedData.gsDistBed] = averageRuns( array2table([df(i).Samples(samplesBed).gsDistMeanFNum], ...
                            'RowNames', df(i).Samples(sampleIdx).gsDistMeanF.Properties.RowNames, ...
                            'VariableNames', matlab.lang.makeUniqueStrings( arrayfun(@(x) x.gsDistMeanF.Properties.VariableNames{:}, ...
                            df(i).Samples(samplesBed), 'Unif', 0)) ) ); % station mean bed composition (unreadable code simply retains row and var names)
                        [df(i).bedData.gsSummBed] = makeSummTable( df(i).bedData.gsDistBed, gsSummQueryPts ); % station mean bed composition
                    end
                catch
                    error('failed to create avg runs, must debug');
                end
            else
                % handler for if there is no bed distribution, just use the
                % mean one from all samples
                df(i).bedData.gsDistBed = meanBedDist;
                df(i).bedData.gsSummBed = meanBedDistSumm;
                df(i).uncertaintyTable.notActualBed = 1;
            end
            % accum water samples
            for j = 1:length(samplesWater)
                sampleIdx = samplesWater(j);
                sampleID = df(i).Samples(sampleIdx).SampleID;
                runIdx = strcmp( sampleIDlist, sampleID );
                if any(runIdx) 
                    % if there is grain size data for this sample specifically, skip to processing
                else
                    warning(sprintf(['no sample distribution found: ', sampleID, '\n', ...
                        'using mean dist of similar sample(s)']))
                    sampleIDStripped = sampleID(1:end-1);
                    runIdx = strcmp( sampleIDlistStripped, sampleIDStripped );
                    df(i).uncertaintyTable.meanWaterDist = 1;
                    if ~any(runIdx)
                        sampleID = 'KT1-C_0.75H_07/03/15a';
                        runIdx = strcmp( sampleIDlist, sampleID );
                        df(i).uncertaintyTable.nearbyWaterDist = 1;
                        df(i).uncertaintyTable.meanWaterDist = 0;
                    end
                end
                if any(runIdx) % if there is a sample distribution found process it
                    [df(i).Samples(sampleIdx)] = accumSampleData(df(i).Samples(sampleIdx), summTable, dataTable, runIdx, gsSummQueryPts);
                    [df(i).Samples(sampleIdx)] = subtractWL(df(i).Samples(sampleIdx), df(i).bedData.gsDistBed, gsPercWLcutoff, gsGSWLcutoff); % grain size distribution without washload
                    [df(i).Samples(sampleIdx).gsSummMeanNWnorm] = makeSummTable(df(i).Samples(sampleIdx).gsDistMeanNWnorm, gsSummQueryPts);
                else
                    x = 1;
                end
            end
            
            % group the samples into table based on near-bed
            samplesNearBed = cell2mat({df(i).Samples.NominalDepth}') == 0.95;
            samplesNearBedDistFCell = arrayfun(@(x) table2array(x.gsDistMeanF), df(i).Samples(samplesNearBed), 'Unif', 0);
            [df(i).nearBedData.gsDistNearBedF] = averageRuns(    array2table(  [samplesNearBedDistFCell{:}], ... % station mean near bed full composition
                    'RowNames', df(i).Samples(samplesNearBed(1)).gsDistMeanF.Properties.RowNames, ...
                    'VariableNames', matlab.lang.makeUniqueStrings( arrayfun(@(x) x.gsDistMeanF.Properties.VariableNames, df(i).Samples(samplesNearBed)) )  )    ); % station mean near bed composition (unreadable code simply retains row and var names)
            samplesNearBedDistNWCell = arrayfun(@(x) table2array(x.gsDistMeanNW), df(i).Samples(samplesNearBed), 'Unif', 0);
            [df(i).nearBedData.gsDistNearBedNW] = averageRuns(    array2table(  [samplesNearBedDistNWCell{:}], ... % station mean near bed composition no wash
                    'RowNames', df(i).Samples(samplesNearBed(1)).gsDistMeanNW.Properties.RowNames, ...
                    'VariableNames', matlab.lang.makeUniqueStrings( arrayfun(@(x) x.gsDistMeanNW.Properties.VariableNames, df(i).Samples(samplesNearBed)) )  )    ); % station mean near bed composition NW (unreadable code simply retains row and var names)
            samplesNearBedDistNWnormCell = arrayfun(@(x) table2array(x.gsDistMeanNWnorm), df(i).Samples(samplesNearBed), 'Unif', 0);
            [df(i).nearBedData.gsDistNearBedNWnorm] = averageRuns(    array2table(  [samplesNearBedDistNWnormCell{:}], ...
                    'RowNames', df(i).Samples(samplesNearBed(1)).gsDistMeanNWnorm.Properties.RowNames, ...
                    'VariableNames', matlab.lang.makeUniqueStrings( arrayfun(@(x) x.gsDistMeanNWnorm.Properties.VariableNames, df(i).Samples(samplesNearBed)) )  )    ); % station mean near bed composition NW (unreadable code simply retains row and var names)
            [df(i).nearBedData.gsSummNearBedNWnorm] = makeSummTable( df(i).nearBedData.gsDistNearBedNWnorm, gsSummQueryPts ); % station mean near bed composition no wash, normalized
            
            % group the samples into water samples table
            [df(i)] = make_waterSamplesTable(df(i));
            
            if simplifyFlag
                % replace simplified versions of datasets
                [df(i)] = make_simpleDatasets(df(i), simplifyNclasses);
            end
            
            if false
                plotStationGrainSize(df(i))
            end
        end
    end
    
    %% save the table
    disp('saving.....')
    save('./dataExport/stationSurveyDataTable.mat', 'df')
    
    
end
