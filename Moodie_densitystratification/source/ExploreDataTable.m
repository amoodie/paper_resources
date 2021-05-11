function ExploreDataTable()
    
    source_path = genpath('source');
    data_path = genpath('datasrc');
    export_path = genpath('dataExport');
    addpath(source_path, data_path, export_path);
    
    load('./dataExport/stationSurveyDataTable.mat', 'df')
    con = load_conset('quartz-water');
    
    WP = readtable('./dataSrc/wright_dissertation_combined.csv');
    WP.Rep = (con.R .* con.g .* (WP.D50_mm_./1e3).^3) ./ con.nu;
    
    if isempty(df(1).velProf)
        error('no velProf listing, are you sure you ran MakeStratModels since running MakeDataTable?')
    end

    if false % recreate the 2018 stations dataset
        df18 = df(arrayfun(@(x) year(x.CollectionDate) == 2018, df), :);
        save('./dataExport/2018_stations.mat', 'df18')
        % where is this used??
    end
    
    global colorSet col1 col2 col3 mar1 mar2 mar3
    colswitch = 'colorful';
    [colorSet, col1, col2, col3, col4] = load_colorSet(colswitch);
    mar1 = 's';
    mar2 = 'o';
    mar3 = '^';
    markers = {mar1; mar2; mar3};
    printOut = false;
    movieOut = false;
    
    mods.WPalphaEvalXs = 1:5:200;
    mods.WPalphaEval = [(1 - (0.06 .* (mods.WPalphaEvalXs(mods.WPalphaEvalXs<=10) .^ 0.77))), ...
        (0.67 - (0.0025 .* (mods.WPalphaEvalXs(mods.WPalphaEvalXs>10))))];
    mods.WPdatacloud = [1.471,0.97293;3.529,0.95726;8.75,0.93875;12.941,0.90456;16.985,0.77778;30.368,0.76211;53.015,0.82621;
                        52.941,0.48148;36.103,0.49858;20.882,0.59259;14.926,0.5641;11.029,0.58262;5.735,0.66952;2.721,0.76923;
                        2.132,0.85897;0.147,0.92593;0,1;1.25,0.97578]; % ripped from plot
    
    [~, nearBed.dims] = tableVarAtDepth(df, 'concNW', 0.95);
    idx.outlierStation = true(size(df, 1), 1);
    idx.outlierStation([18 24 47 31 49 12 52]) = false;
%     outlierStation = arrayfun(@(x) ismember(x.WrtChannel, {'C', 'LB'}), df); % all but center are outliers
%     idx.sample = logical(  varRepLong( double(arrayfun(@(x) strcmp(x.NominalLocation, 'Kenli'), df)), nearBedDims )  ); % highlighting index for samples
    idx.sampleYr = varRepLong( double(arrayfun(@(x) year(datenum(x.CollectionDate)), df)), nearBed.dims ); % highlighting index for samples
    idx.stationYr = double(arrayfun(@(x) year(datenum(x.CollectionDate)), df)); % highlighting index for stations
    idx.stationOut = and(idx.stationYr, idx.outlierStation);
    idx.sampleOut = varRepLong(double(idx.stationOut), nearBed.dims);
    
    velProfIdx = arrayfun(@(x) ~(size(fieldnames(x.velProf),1)<5), df);
    nClasses = size(df(1).nearBedData.gsDistNearBedNWnorm, 1);
    idx.velProf = arrayfun(@(x) ~(size(fieldnames(x.velProf),1)<5), df);
    
    %% nearBed.concNW vs (ustar, depth, ustar/ws*Rep) all data
    [nearBed.concWash, ~] = tableVarAtDepth(df, 'conc', 0.95); % varAtDepth(df, varName, depth)
    [nearBed.concNW, ~] = tableVarAtDepth(df, 'concNW', 0.95);
    [nearBed.discharge] = varRepLong(arrayfun(@(x) x.Velocity.discharge, df), nearBed.dims);
    [nearBed.concNWVol] = nearBed.concNW ./ 2650;
    [nearBed.concNWbyClass, ~] = tableVarAtDepth(df, 'concNWbyClass', 0.95);
    [nearBed.concNWVolbyClass] = concatenate_cell_vectors(nearBed.concNWbyClass) ./ 2650;
    [nearBed.NWd50] = tableVarAtDepth(df, 'gsSummNWnorm.d50', 0.95);
    [nearBed.NWd10] = tableVarAtDepth(df, 'gsSummNWnorm.d10', 0.95);
    [nearBed.NWd90] = tableVarAtDepth(df, 'gsSummNWnorm.d90', 0.95);
    [nearBed.NWd50ws] = get_DSV(nearBed.NWd50 ./ 1e6, 0.7, 3.5, con);
    [nearBed.FlowDepth] = varVsDepth(df, 'FlowDepthAtCollection', nearBed.dims);
    [nearBed.ustarDSP] = varVsDepth(df, 'Velocity.ustarDSP', nearBed.dims);
    [nearBed.ustarCalib] = varVsDepth(df, 'Velocity.ustarCalib', nearBed.dims);
    [nearBed.ustarMeas] = varVsDepth(df, 'Velocity.ustarCalib', nearBed.dims);
    [nearBed.ustarsk] = varVsDepth(df, 'Velocity.ustarsk', nearBed.dims);
    [nearBed.Ubar_conCf] = varVsDepth(df, 'Velocity.Ubar_conCf', nearBed.dims);
    [nearBed.Rep] = ( sqrt(con.R * con.g .* nearBed.NWd50 .* 1e-6) .* nearBed.NWd50.*1e-6 ) ./ con.nu;
    
    [bed.d10] = varVsDepth(df, 'bedData.gsSummBed.d10', nearBed.dims);
    [bed.d50] = varVsDepth(df, 'bedData.gsSummBed.d50', nearBed.dims);
    [bed.d90] = varVsDepth(df, 'bedData.gsSummBed.d90', nearBed.dims);
    [bed.d10ws] = get_DSV(bed.d10 ./ 1e6, 0.7, 3.5, con);
    [bed.d50ws] = get_DSV(bed.d50 ./ 1e6, 0.7, 3.5, con);
    [bed.d90ws] = get_DSV(bed.d90 ./ 1e6, 0.7, 3.5, con);
    [bed.d50Rep] = (sqrt(con.R .* con.g .* (bed.d50*1e-6)) .* (bed.d50*1e-6)) ./ con.nu;
    [bed.slope] = varVsDepth(df, 'Velocity.slope', nearBed.dims);
    
%     [d50ws] = bed.d50ws; % this can be switched to near bed for testing (?)
    
    station.year = arrayfun(@(x) year(x.CollectionDate), df);
    
    station.tstar = arrayfun(@(x) x.Velocity.tstar, df); 
    station.ustarCalib = arrayfun(@(x) x.Velocity.ustarCalib, df); 
    station.ustarDSP = arrayfun(@(x) x.Velocity.ustarDSP, df);
    station.ustarsk = arrayfun(@(x) x.Velocity.ustarsk, df);
    station.depth = arrayfun(@(x) x.FlowDepthAtCollection, df);
    station.slope = arrayfun(@(x) x.Velocity.slope, df); % determined slope
    station.Ubar = arrayfun(@(x) (mean(x.velProf.MYfullPred.Us)), df);
    station.Cm = arrayfun(@(x) x.r0calcs.Cm, df); % discharge weighted total concentration by volume        
    station.Fr = arrayfun(@(x) mean(x.velProf.MYfullPred.Us)./sqrt(con.g * x.FlowDepthAtCollection), df); % discharge weighted total concentration by volume
    station.Cf = arrayfun(@(x) (con.g * x.FlowDepthAtCollection * x.Velocity.slope) / (mean(x.velProf.MYfullPred.Us)^2), df);
    
    station.suspension.nearBedNWnormD50ws = arrayfun(@(x) get_DSV(x.nearBedData.gsSummNearBedNWnorm.d50 ./ 1e6, 0.7, 3.5, con), df); 
    station.suspension.nearBedNWnormDist = arrayfun(@(x) gsDistClass2num(x.nearBedData.gsDistNearBedNWnorm), df, 'Unif', 0);
    station.suspension.nearBedNWnormDistVals = arrayfun(@(x) table2array(x.nearBedData.gsDistNearBedNWnorm), df, 'Unif', 0);
    station.suspension.nearBedNWnormDistWs = cellfun(@(x) get_DSV(x ./ 1e6, 0.7, 3.5, con), station.suspension.nearBedNWnormDist, 'Unif', 0); 
    station.suspension.nearBedNWnormDistWeightedWs = cellfun(@(dist, val) nansum((dist ./ 100) .* val), station.suspension.nearBedNWnormDistVals, station.suspension.nearBedNWnormDistWs);
    station.suspension.suspensionNumberNearBedNWnormD50 = station.ustarCalib ./ station.suspension.nearBedNWnormD50ws;
    station.suspension.suspensionNumberNearBedNWnormDistWeighted = station.ustarCalib ./ station.suspension.nearBedNWnormDistWeightedWs;
    station.suspension.nearBedNWnormDistRep = cellfun(@(x) ( sqrt(con.R * con.g .* (x./1e6)) .* (x./1e6) ./ con.nu), station.suspension.nearBedNWnormDist, 'Unif', 0); 
    station.suspension.nearBedNWnormRepD50 = ( sqrt(con.R * con.g .* arrayfun(@(x) x.nearBedData.gsSummNearBedNWnorm.d50, df).*1e-6) ...
        .* arrayfun(@(x) x.nearBedData.gsSummNearBedNWnorm.d50, df).*1e-6 ) ./ con.nu;
    station.suspension.nearBedNWnormDistWeightedRep = cellfun(@(dist, val) nansum((dist ./ 100) .* val), arrayfun(@(x) table2array(x.nearBedData.gsDistNearBedNWnorm), df, 'Unif', 0), station.suspension.nearBedNWnormDistRep);
    
    station.suspension.bedD50 = arrayfun(@(x) x.bedData.gsSummBed.d50, df) .* 1e-6;
    station.suspension.bedD50star = ((con.R * con.g) / (con.nu * con.nu))^(1/3) .* station.suspension.bedD50;
    station.suspension.bedD50ws = arrayfun(@(x) get_DSV(x.bedData.gsSummBed.d50 ./ 1e6, 0.7, 3.5, con), df); 
    station.suspension.bedDist = arrayfun(@(x) gsDistClass2num(x.bedData.gsDistBed), df, 'Unif', 0);
    station.suspension.bedDistWs = cellfun(@(x) get_DSV(x ./ 1e6, 0.7, 3.5, con), station.suspension.bedDist, 'Unif', 0); 
    station.suspension.bedDistWeightedWs = cellfun(@(dist, val) nansum((dist ./ 100) .* val), ...
        arrayfun(@(x) table2array(x.bedData.gsDistBed), df, 'Unif', 0), station.suspension.bedDistWs); 
    station.suspension.suspensionNumberBedD50 = station.ustarCalib ./ station.suspension.bedD50ws;
    station.suspension.suspensionNumberBedDistWeighted = station.ustarCalib ./ station.suspension.bedDistWeightedWs;
    station.suspension.bedDistRep = cellfun(@(x) ( sqrt(con.R * con.g .* (x./1e6)) .* (x./1e6) ./ con.nu), station.suspension.bedDist, 'Unif', 0); 
    station.suspension.bedRepD50 = ( sqrt(con.R * con.g .* arrayfun(@(x) x.bedData.gsSummBed.d50, df).*1e-6) ...
        .* arrayfun(@(x) x.bedData.gsSummBed.d50, df).*1e-6 ) ./ con.nu;
    station.suspension.bedDistWeightedRep = cellfun(@(dist, val) nansum((dist ./ 100) .* val), arrayfun(@(x) table2array(x.bedData.gsDistBed), df, 'Unif', 0), station.suspension.bedDistRep);    
    
    station.Rou.gsClassPredD50Rou = arrayfun(@(x) x.concProf.gsClassPred.d50Rou, df);
    station.Rou.totalModelRou = arrayfun(@(x) x.concProf.totalModel.params.Rou, df);
    station.Rou.gsClassPredBulkRou = arrayfun(@(x) x.concProf.gsClassPred.bulkRou, df);
    station.Rou.gsClassModelBulkRou = arrayfun(@(x) x.concProf.gsClassModel.bulkRou, df);
    
    [station.theta.byDepthToGsClassModelNorm] = arrayfun(@(x) x.concProf.theta.byDepthToGsClassModel, df, 'unif', 0);
        [station.theta.byDepthToGsClassModelNorm] = [station.theta.byDepthToGsClassModelNorm{:}]';
    [station.theta.byDepthMYfullGsClassToGsClassModelNorm] = arrayfun(@(x) nansum(x.concProf.theta.byDepthMYfullGsClassToGsClassModel, 2), df, 'unif', 0);
        [station.theta.byDepthMYfullGsClassToGsClassModelNorm] = [ station.theta.byDepthMYfullGsClassToGsClassModelNorm{:} ]';
    [station.theta.byDepthMYfullGsClassToGsClassModelDepthNorm] = arrayfun(@(x) x.concProf.theta.byDepthMYfullGsClassToGsClassModelNorm, df, 'unif', 0);
        [station.theta.byDepthMYfullGsClassToGsClassModelDepthNorm] = cat(3, station.theta.byDepthMYfullGsClassToGsClassModelDepthNorm{:});
    [station.theta.meanNormToTotalModel] = arrayfun(@(x) x.concProf.theta.meanNormToTotalModel, df);
    [station.theta.meanNormToGsClassModel] = arrayfun(@(x) x.concProf.theta.meanNormToGsClassModel, df);
    [station.theta.meanToTotalModel] = arrayfun(@(x) x.concProf.theta.meanToTotalModel, df);
    [station.theta.meanToGsClassModel] = arrayfun(@(x) x.concProf.theta.meanToGsClassModel, df);
    
    
%     [station.nearBed.concNW] = arrayfun(@(x) x.concProf.gsClassPred.CsSum(1), df);
%     
%     station.bed.Rep
%     [station.bed.Repd50] = ( sqrt(con.R * con.g .* arrayfun(@(x) x.bedData.gsSummBed.d50, df).*1e-6) ...
%         .* arrayfun(@(x) x.bedData.gsSummBed.d50, df).*1e-6 ) ./ con.nu;
%     [station.bed.RepBulk] = cellfun(@(dist, val) nansum((dist ./ 100) .* val), ...
%         arrayfun(@(x) table2array(x.bedData.gsDistBed), df, 'Unif', 0), station.suspension.bedDistWs);
      
    
    station.Risr.byBedD50 = con.R .* (station.suspension.bedD50ws ./ station.ustarCalib) .* (station.Cm ./ station.slope); % sand richardson number
    station.Risr.byBedBulk = con.R .* (station.suspension.bedDistWeightedWs ./ station.ustarCalib) .* (station.Cm ./ station.slope); % sand richardson number
    station.Risr.byNearBedBulk = con.R .* (station.suspension.nearBedNWnormDistWeightedWs ./ station.ustarCalib) .* (station.Cm ./ station.slope); % sand richardson number
    
    [station.alpha.byPredBulkRouToTotalModelRou] = double(arrayfun(@(x) x.concProf.alpha.byPredBulkRouToTotalModelRou, df));
    [station.alpha.byPredBulkRouToGsClassModelRou] = concatenate_cell_vectors( arrayfun(@(x) x.concProf.alpha.byPredBulkRouToGsClassModelRou, df, 'unif',0) );
    [station.alpha.byPredD50RouToTotalModelRou] = double(arrayfun(@(x) x.concProf.alpha.byPredD50RouToTotalModelRou, df));
    [station.alpha.byPredD50RouToGsClassModelRou] = concatenate_cell_vectors( arrayfun(@(x) x.concProf.alpha.byPredD50RouToGsClassModelRou, df, 'unif',0) );
    [station.alpha.byResid] = double(arrayfun(@(x) x.concProf.alpha.byResid, df));
    [station.alpha.byKredMY] = double(arrayfun(@(x) x.concProf.alpha.byKredMY, df));
    [station.alpha.byWP04func] = double(arrayfun(@(x) x.concProf.alpha.byWP04func, df));
    [station.alpha.byKredVelocityModel] = arrayfun(@(x) x.velProf.modelKred.params.alpha, df(velProfIdx));
    
    [station.cb5S.byTotalModel] = ( arrayfun(@(x) x.concProf.totalModel.Cs(1), df) ./ 2.65 ./ 1000) ./ arrayfun(@(x) x.Velocity.slope, df); % cb/S in Wright and Parker alpha relation
    [station.cb5S.byGsClassSummedModel] = ( arrayfun(@(x) x.concProf.gsClassModel.CsSum(1), df) ./ 2.65 ./ 1000) ./ arrayfun(@(x) x.Velocity.slope, df); % cb/S in Wright and Parker alpha relation
    [station.cb5S.cb5] = (arrayfun(@(x) x.concProf.gsClassModel.CsSum(1), df) ./ 2.65 ./ 1000);
    [station.cb5S.S] = arrayfun(@(x) x.Velocity.slope, df);
     
    % grain size specific station lists
    [station.theta.meanNormGsClassToGsClassModel] = concatenate_cell_vectors( arrayfun(@(x) vertcat(x.concProf.theta.meanNormGsClassToGsClassModel), df, 'unif', 0) ); 
    [station.theta.meanGsClassToGsClassModel] = concatenate_cell_vectors( arrayfun(@(x) vertcat(x.concProf.theta.meanGsClassToGsClassModel), df, 'unif', 0) );
    [station.alpha.byPredGsClassRouToGsClassModelRou] = concatenate_cell_vectors( arrayfun(@(x) x.concProf.alpha.byPredGsClassRouToGsClassModelRou, df, 'unif',0) );
    [station.gsClassVector] = concatenate_cell_vectors( arrayfun(@(x) vertcat(x.concProf.gsClassModel.gsClass), df, 'unif', 0) );
    % [station.gsClassDistAboveBed] = concatenate_cell_vectors( arrayfun(@(x) x.waterSamplesTable.sampleZnorm, df, 'unif', 0) );
    [station.gsClassVectorCats] = categorical(cellstr(num2str(round(station.gsClassVector))));
    [station.cb5S.byGsClassModel] = concatenate_cell_vectors( arrayfun(@(x) x.concProf.gsClassModel.Cs(1,:) ./ 2.65 ./ 1000, df, 'unif', 0) ) ./ ...
        varRepLong(arrayfun(@(x) x.Velocity.slope, df), length(df(1).concProf.gsClassModel.gsClass)); % cb/S in Wright and Parker alpha relation
    [station.cb5S.byGsClassModelToGsClassSummed] = varRepLong(station.cb5S.byGsClassSummedModel, length(df(1).concProf.gsClassModel.gsClass)); % cb/S in Wright and Parker alpha relation
    [station.cb5S.byGsClassModelToGsClassFrac] = station.cb5S.byGsClassModel ./ (concatenate_cell_vectors(station.suspension.nearBedNWnormDistVals)./100) ;
    [station.suspension.gsClassWs] = get_DSV(station.gsClassVector ./ 1e6, 0.7, 3.5, con); 
    [station.suspension.suspensionNumberGsClass] = varRepLong(station.ustarCalib,  length(df(1).concProf.gsClassModel.gsClass)) ./ station.suspension.gsClassWs;
    [station.transport.field_Rou] = concatenate_cell_vectors(arrayfun(@(x) x.concProf.transport.field_Rou, df, 'unif', 0));     
    [station.transport.pred_Rou] = concatenate_cell_vectors(arrayfun(@(x) x.concProf.transport.pred_Rou, df, 'unif', 0)); 
    [station.suspension.nearBedtstarGsClass] = concatenate_cell_vectors( arrayfun(@(x) (x.Velocity.ustarCalib .* x.Velocity.ustarCalib) ./ (9.81 .* 1.65 .* x.concProf.gsClassModel.gsClass.*1e-6), df, 'unif', 0) );
    [station.suspension.suspensionNumberUstarskGsClass] = varRepLong(station.ustarsk,  length(df(1).concProf.gsClassModel.gsClass)) ./ station.suspension.gsClassWs;
    [station.transport.HMGEHGsClass] = sqrt(con.R.*con.g.*(station.gsClassVector.*1e-6).^3)  .* ((0.90 .* station.suspension.nearBedtstarGsClass.^(1.68))./varRepLong(station.Cf, nClasses)) .* (concatenate_cell_vectors(station.suspension.nearBedNWnormDistVals)./100);
    AiKN = @(DiDg) 0.455.*(DiDg).^(-0.839); BiKN = @(DiDg) 0.353.*(DiDg).^(-1.157);
    [station.transport.DiDg] = station.gsClassVector.*1e-6 ./ repelem(station.suspension.bedD50, nClasses);
    [station.transport.KNGsClass] = (concatenate_cell_vectors(station.suspension.nearBedNWnormDistVals)./100) .* ...
        (repelem(station.ustarCalib.^3, nClasses)./(repelem(station.Cf, nClasses).*con.R.*con.g)) .* ...
        ((AiKN(station.transport.DiDg) .* station.suspension.nearBedtstarGsClass .^ (BiKN(station.transport.DiDg))));
    [station.transport.qsMYmodel] = concatenate_cell_vectors(arrayfun(@(x) x.concProf.transport.MYmodelGsClass, df, 'unif', 0));
    
    % RMSE metrics
    station.RMSE.gsClassModel = arrayfun(@(x) x.concProf.gsClassModel.RMSE, df);
    station.RMSE.gsClassModelnorm = arrayfun(@(x) x.concProf.gsClassModel.RMSEnorm, df);
    station.RMSE.gsClassPred = arrayfun(@(x) x.concProf.gsClassPred.RMSE, df);
    station.RMSE.gsClassPrednorm = arrayfun(@(x) x.concProf.gsClassPred.RMSEnorm, df);
    station.RMSE.WPalphaPred = arrayfun(@(x) x.concProf.WPalphaPred.RMSE, df);
    station.RMSE.WPalphaPrednorm = arrayfun(@(x) x.concProf.WPalphaPred.RMSEnorm, df);
    station.RMSE.MYfullPred = arrayfun(@(x) x.concProf.MYfullPred.RMSE, df);
    station.RMSE.MYfullPred_class = concatenate_cell_vectors(arrayfun(@(x) x.concProf.MYfullPred.class_RMSE, df, 'Unif', 0));
    station.RMSE.MYfullPrednorm = arrayfun(@(x) x.concProf.MYfullPred.RMSEnorm, df);
    station.RMSE.MYfullPrednorm_class = concatenate_cell_vectors(arrayfun(@(x) x.concProf.MYfullPred.class_RMSEnorm, df, 'Unif', 0));
    station.RMSE.dataTable = array2table([station.RMSE.MYfullPrednorm, station.RMSE.WPalphaPrednorm, station.RMSE.gsClassPrednorm], 'VariableNames', ...
                                         {'MYfullPredNorm', 'WPalphaPrednorm' , 'gsClassPrednorm'});
    station.RMSE.summaryTable = array2table([nanmean(table2array(station.RMSE.dataTable), 1);  std(table2array(station.RMSE.dataTable), 1, 'omitnan')], ...
                                            'VariableNames', station.RMSE.dataTable.Properties.VariableNames, 'RowNames', {'mean', 'std'});
    [H, P] = ttest2(station.RMSE.dataTable.MYfullPredNorm, station.RMSE.dataTable.WPalphaPrednorm);
    
    station.beta.long.cellmat = arrayfun(@(x) x.concProf.MYfullPred.class_E, df, 'Unif', 0);
    station.beta.long.E = concatenate_cell_vectors( cellfun(@(x) reshape(transpose(x),[numel(x), 1]), station.beta.long.cellmat, 'Unif', 0)  );
    station.beta.long.GsClass = concatenate_cell_vectors( cellfun(@(x) reshape(transpose(x),[numel(x), 1]), arrayfun(@(x) x.concProf.MYfullPred.class_class, df, 'Unif', 0), 'Unif', 0)  );
    station.beta.long.DistNorm = concatenate_cell_vectors( cellfun(@(x) reshape(transpose(x),[numel(x), 1]), arrayfun(@(x) x.concProf.MYfullPred.class_DistNorm, df, 'Unif', 0), 'Unif', 0)  );
    
    % r0 and sed transport calcs
    station.beta.fieldr0 = arrayfun(@(x) x.concProf.beta.fieldr0, df, 'Unif', 0);
    station.beta.Predr0 = arrayfun(@(x) x.concProf.beta.Predr0, df, 'Unif', 0);
    station.beta.fieldr0Long = concatenate_cell_vectors(station.beta.fieldr0);
    station.beta.Predr0Long = concatenate_cell_vectors(station.beta.Predr0);
    station.beta.fieldr0Long(station.beta.fieldr0Long < 1e-3) = NaN; station.beta.Predr0Long(station.beta.Predr0Long < 1e-3) = NaN;
    station.beta.alphaGsClassVelProfOnly = concatenate_cell_vectors( arrayfun(@(x) x.concProf.alpha.byPredGsClassRouToGsClassModelRou, df(idx.velProf), 'unif',0) );
    station.beta.nearBedConcVelProfOnly = concatenate_cell_vectors(arrayfun(@(x) x.concProf.gsClassModel.Cs(1,:), df(idx.velProf), 'unif', 0));
    station.beta.alphaByVelocityVelProfOnly = repelem(station.alpha.byKredVelocityModel, nClasses);
    station.beta.byGsClassModelToGsClassFrac = station.beta.nearBedConcVelProfOnly ./ (concatenate_cell_vectors(station.suspension.nearBedNWnormDistVals(idx.velProf))./100) ;
    
    station.transport.MYmodel = arrayfun(@(x) x.concProf.transport.MYmodel, df);
    station.transport.noStrat = arrayfun(@(x) x.concProf.transport.noStrat, df);
    station.transport.noStrat_upper = arrayfun(@(x) x.concProf.transport.noStrat_upper, df);
    station.transport.field = arrayfun(@(x) x.concProf.transport.fits, df(velProfIdx));
    station.transport.field_upper = arrayfun(@(x) x.concProf.transport.fits_upper, df(velProfIdx));
%     station.transport.fieldAll = arrayfun(@(x) x.concProf.transport.fits, df);
%     station.transport.fieldr0Bulk = arrayfun(@(x) x.concProf.beta.fieldr0Bulk, df);
    station.transport.fieldr0Bulk = arrayfun(@(x) x.concProf.gsClassModel.CsSum(1,:) ./ mean(x.concProf.gsClassModel.CsSum, 1), df);
    station.transport.MYnoStrat_ratio = station.transport.MYmodel ./ station.transport.noStrat;
    station.transport.fieldNoStrat_ratio = station.transport.field ./ station.transport.noStrat(velProfIdx);
    station.transport.fieldNoStrat_ratio_upper = station.transport.field_upper ./ station.transport.noStrat_upper(velProfIdx);
    station.transport.pred_D10 = arrayfun(@(x) x.concProf.transport.pred_D_ii(1), df);
    station.transport.fit_D10 = arrayfun(@(x) x.concProf.transport.fit_D_ii(1), df);
    station.transport.pred_D50 = arrayfun(@(x) x.concProf.transport.pred_D_ii(2), df);
    station.transport.fit_D50 = arrayfun(@(x) x.concProf.transport.fit_D_ii(2), df);
    station.transport.pred_D90 = arrayfun(@(x) x.concProf.transport.pred_D_ii(3), df);
    station.transport.fit_D90 = arrayfun(@(x) x.concProf.transport.fit_D_ii(3), df);
    station.transport.pred_sandvol = arrayfun(@(x) nansum(x.concProf.transport.pred_perconcbin(2:end)), df);
    station.transport.fit_sandvol = arrayfun(@(x) nansum(x.concProf.transport.fit_perconcbin(2:end)), df);
    station.transport.concRatio = arrayfun(@(x) x.concProf.transport.concRatio, df);
    station.transport.concRatioGsClass = concatenate_cell_vectors(arrayfun(@(x) x.concProf.transport.concRatioGsClass, df, 'unif', 0));
    station.transport.concFracReduct = arrayfun(@(x) x.concProf.transport.concFracReduct, df);
    station.transport.concFracReductSumm = [mean(station.transport.concFracReduct(idx.stationOut)), std(station.transport.concFracReduct(idx.stationOut))];
    [station.transport.HMGEH] = sqrt(con.R.*con.g.*(station.suspension.bedD50).^3) .* ((0.90 .* station.tstar.^(1.68)) ./ station.Cf);

    [station.transport.Esolve] = ((station.transport.HMGEHGsClass .* station.beta.fieldr0Long)) ./ repelem(station.depth .* station.Ubar, nClasses);
    
    
    % velocity calibration
    station.velCalib.depth = arrayfun(@(x) x.FlowDepthAtCollection, df(velProfIdx));
    station.velCalib.slope = arrayfun(@(x) x.Velocity.slope, df(velProfIdx));
    station.velCalib.discharge = arrayfun(@(x) x.Velocity.discharge, df(velProfIdx));
    station.velCalib.ustarsk = arrayfun(@(x) x.Velocity.ustarsk, df(velProfIdx));
    station.velCalib.Ubar_fromdata = arrayfun(@(x) x.Velocity.data.mean, df(velProfIdx));
    station.velCalib.ustarFric = arrayfun(@(x) x.Velocity.ustarFric, df(velProfIdx));
    station.velCalib.ustarFric_MaSolved = arrayfun(@(x) x.Velocity.ustarFric_MaSolved, df(velProfIdx));
    station.velCalib.ustarDSP = arrayfun(@(x) x.Velocity.ustarDSP, df(velProfIdx));
    station.velCalib.ustarNoAlpha = arrayfun(@(x) x.velProf.modelNoAlpha.params.ustar, df(velProfIdx));
    station.velCalib.ustarAlpha = arrayfun(@(x) x.velProf.modelAlpha.params.ustar, df(velProfIdx));
    station.velCalib.ustarYr = double(arrayfun(@(x) year(datenum(x.CollectionDate)), df(velProfIdx)));   
    velocityCalibTable = array2table([[station.velCalib.discharge], ...
                                 [station.velCalib.depth], ...
                                 [station.velCalib.slope], ...
                                 [station.velCalib.ustarDSP], ...
                                 [station.velCalib.ustarFric], ...
                                 [station.velCalib.ustarNoAlpha], ...
                                 [station.velCalib.ustarAlpha], ...
                                 [station.velCalib.ustarYr]], 'VariableNames', ...
                                 {'discharge', 'depth', 'slope', 'ustarDSP', 'ustarFric', ...
                                 'ustarNoAlpha', 'ustarAlpha', 'year'});
                             
                             
    % reduce by outliers and make data ready for regression analysis
    dfd = df(idx.outlierStation); % drop the outliers
    wss = station.suspension.nearBedNWnormDistWs(idx.outlierStation);
    reg.concs = NaN([size(dfd(1).concProf.gsClassModel.Cs), size(dfd,1)]);
    reg.ustarwss = NaN([1, size(dfd(1).concProf.gsClassModel.Cs, 2), size(dfd,1)]);
    reg.cbiSs = NaN([1, size(dfd(1).concProf.gsClassModel.Cs, 2), size(dfd,1)]);
    reg.cbis = NaN([1, size(dfd(1).concProf.gsClassModel.Cs, 2), size(dfd,1)]);
    reg.cbSs = NaN([1, 1, size(dfd,1)]);
    reg.S = NaN([1, 1, size(dfd,1)]);
    for i = 1:size(dfd,1)
        reg.concs(:,:,i) = dfd(i).concProf.gsClassModel.Cs ./ 2650;
        reg.ustarwss(1,:,i) = dfd(i).Velocity.ustarCalib ./ wss{i}';
        reg.cbis(1,:,i) = (dfd(i).concProf.gsClassModel.Cs(1,:) ./ 2650);
        reg.cbiSs(1,:,i) = (dfd(i).concProf.gsClassModel.Cs(1,:) ./ 2650) ./ dfd(i).Velocity.slope;
        reg.cbSs(1, 1, i) = (dfd(i).concProf.gsClassModel.CsSum(1) ./ 2650) ./ dfd(i).Velocity.slope;
        reg.S(1, 1, i) = dfd(i).Velocity.slope;
    end
    reg.Zs = linspace(0.05, 1, 51)';
    reg.gamma0 = ones(size(reg.ustarwss));
    reg.alphas0 = reg.gamma0 .* calculate_alpha_WP04(sum(reg.cbis, 2), reg.S);
    reg.rouses0 = (1./reg.ustarwss) .* (1./(0.41 .* reg.alphas0));
    reg.pred_concs0 = reg.cbis .* ( ((1-reg.Zs)./reg.Zs).*(0.05./(1-0.05)) ) .^ reg.rouses0;
    reg.msd = mean(reg.pred_concs0(3:end-2,:,:) - reg.concs(3:end-2,:,:), 1); % mean signed deviation for every profile, when normalized in 0-1 cb space
    reg.nmsd = reg.msd ./ reg.cbis; % normalized to nearbed
    reg.sse = ( sum( (reg.pred_concs0-reg.concs).^2, 1) / size(reg.pred_concs0,1) ).^(0.5);
    reg.sdd = nanmean(reg.pred_concs0 - reg.concs, 3); %./  mean(cbis,3);
    reg.sddstd = std(reg.pred_concs0 - reg.concs, [], 3, 'omitnan');
    
    % do the regression for correction Gamma
%     reg.corr = 
    
    
    % set up entrainment stuff
    %% entrainment 'parker' relations
    entr.mods.GP91data = array2table(csvread('./dataSrc/GP1991_data.csv', 0, 0), 'VariableNames', {'Zu', 'Es'});
    entr.mods.GP91Zu = logspace(0, 2, 100); 
    entr.mods.GP91A = 1.3 * 10^-7;
    entr.mods.GP91Es = (entr.mods.GP91A .* entr.mods.GP91Zu .^ 5) ./ (1 + ((entr.mods.GP91A / 0.3) .* entr.mods.GP91Zu .^ 5));
    entr.mods.WP04Zu = entr.mods.GP91Zu;
    entr.mods.WP04B = 7.8 * 10^-7;
    entr.mods.WP04Es = (entr.mods.WP04B .* entr.mods.WP04Zu .^ 5) ./ (1 + ((entr.mods.WP04B / 0.3) .* entr.mods.WP04Zu .^ 5));
    entr.mods.WP04EsNoLim = (entr.mods.WP04B .* entr.mods.WP04Zu .^ 5);
        
    entr.d50.GP91_Zu = (nearBed.ustarDSP ./ bed.d50ws) .* (bed.d50Rep .^ 0.6);
    entr.d50.WP04_Zu = (nearBed.ustarDSP ./ bed.d50ws) .* (bed.d50Rep .^ 0.6) .* (bed.slope .^ 0.08);
   
    Slope = 6.4e-5;
    slope_exp = 0.08;
    entr.Slope = Slope;
    entr.slope_exp = slope_exp;
    
    [entr.dist.cellAll] = arrayfun(@(s) {eval(['s.' 'bedData.gsDistBed'])}, df, 'Unif', 0);
    [entr.dist.cellClass] = cellfun(@(x) str2double(x{1}.Properties.RowNames'), entr.dist.cellAll, 'Unif', 0);
    [entr.dist.cellDist] = cellfun(@(x) x{1}.avgClassNorm', entr.dist.cellAll, 'Unif', 0);
    [entr.dist.dist] = varRepLong(entr.dist.cellDist, nearBed.dims);
    [entr.dist.class] = varRepLong(entr.dist.cellClass, nearBed.dims);
    [entr.dist.ustarsk] = varRepLong(vertcat(arrayfun(@(s) eval(['s.' 'Velocity.ustarCalib']), df, 'Unif', 0)), nearBed.dims); % varRepLong(vertcat(arrayfun(@(s) eval(['s.' 'Velocity.ustarCalib']), df, 'Unif', 0)), nearBed.dims);
    [entr.dist.ws] = cellfun(@(class) get_DSV(class.*1e-6, 0.7, 3.5, con), entr.dist.class, 'Unif', 0);
    [entr.dist.Rep] = cellfun(@(class) (sqrt(con.R .* con.g .* (class.*1e-6)) .* (class.*1e-6)) ./ con.nu, entr.dist.class, 'Unif', 0);
    [entr.dist.ustarskWs] = cellfun(@(ustarsk, ws) ustarsk ./ ws, entr.dist.ustarsk, entr.dist.ws, 'Unif', 0);
    [entr.dist.DiD50] = cellfun(@(D50, class) (class./D50), num2cell(bed.d50), entr.dist.class, 'Unif', 0); % Di/D50
    [entr.dist.ustarskWsRep1] = cellfun(@(ustarsk, ws, Rep) ustarsk ./ ws .* Rep .^ 1, entr.dist.ustarsk, entr.dist.ws, entr.dist.Rep, 'Unif', 0);
    [entr.dist.lam] = varRepLong(arrayfun(@(s) eval(['s.' 'bedData.gsSummBed.lambda']), df), nearBed.dims);
    [entr.dist.Xi] = cellfun(@(ustarsk, ws, Rep, D50, class) (((ustarsk ./ ws).*(Rep.^0.6)).*Slope.^(slope_exp).*(((class./D50)).^0.2)), ...
        (entr.dist.ustarsk), entr.dist.ws, entr.dist.Rep, num2cell(bed.d50), entr.dist.class, 'Unif', 0); % ustar sk
    [entr.dist.ustarskWsRep06] = cellfun(@(ustarsk, ws, Rep) ustarsk ./ ws .* Rep .^ 0.6, entr.dist.ustarsk, entr.dist.ws, entr.dist.Rep, 'Unif', 0);
    [entr.dist.ustarskWSRep06long] = concatenate_cell_vectors(entr.dist.ustarskWsRep06);
    [entr.dist.XiEff] = cellfun(@(Xi, d) trapz(cumsum(d./100), Xi), entr.dist.Xi, entr.dist.dist);
    [entr.dist.XiWeight] = cellfun(@(Xi, d) nansum(Xi .* (d./100)), entr.dist.Xi, entr.dist.dist);
    [entr.dist.lXi] = cellfun(@(lam, Xi) lam .* Xi, num2cell(entr.dist.lam), entr.dist.Xi, 'Unif', 0);
    [entr.dist.lXiLong] = concatenate_cell_vectors(entr.dist.lXi);
    %[entr.dist.lXiSum] = cellfun(@(dist, lam, Xi) lam * nansum((dist ./ 100) .* Xi), entr.dist.dist, num2cell(entr.dist.lam), entr.dist.Xi);
    [entr.dist.lXiEff] = cellfun(@(lam, XiEff) lam .* XiEff, num2cell(entr.dist.lam), num2cell(entr.dist.XiEff));
    [entr.dist.lXiWeight] = cellfun(@(lam, XiWeight) lam .* XiWeight, num2cell(entr.dist.lam), num2cell(entr.dist.XiWeight));
    [entr.dist.EsiPred] = cellfun(@(lXi, Fi) ( (7.8e-7 .* (lXi) .^ 5) ./ (1 + (7.8e-7/0.3) .* (lXi).^5) ) , ...
        entr.dist.lXi, entr.dist.dist, 'Unif', 0);
    [entr.dist.EsiMeas] = cellfun(@(data, Fi) (data./2650)' ./ (Fi./100), nearBed.concNWbyClass, entr.dist.dist, 'Unif', 0);
    [entr.dist.EsiMeasNotNorm] = cellfun(@(data) (data./2650), nearBed.concNWbyClass, 'Unif', 0);
%     [entr.dist.Esilong] = concatenate_cell_vectors(entr.dist.Esi);
    [entr.dist.EsPred] = cellfun(@(Esi, Fi) nansum(Esi(1:end) .* (Fi(1:end)./100)), entr.dist.EsiPred, entr.dist.dist);
    [entr.dist.EsMeas] = nearBed.concNW./2650; % cellfun(@(Esi, Fi) nansum(Esi(~isinf(Esi)) .* (Fi(~isinf(Esi))./100)), entr.dist.EsiMeas, entr.dist.dist); 
    [entr.dist.EsiPredLong] = concatenate_cell_vectors(entr.dist.EsiPred);
    [entr.dist.EsiMeasLong] = concatenate_cell_vectors(entr.dist.EsiMeas);
    [entr.dist.EsiMeasNotNormLong] = concatenate_cell_vectors(entr.dist.EsiMeasNotNorm);
    [entr.dist.EsiPredLongCourseOnly] = concatenate_cell_vectors(cellfun(@(Esi) Esi(2:end), entr.dist.EsiPred, 'Unif', 0));
    [entr.dist.EsiMeasLongCourseOnly] = concatenate_cell_vectors(cellfun(@(Esi) Esi(2:end), entr.dist.EsiMeas, 'Unif', 0));
                    
    %% sediment diffusivity metrics
    % reduce by outliers and make data ready for regression analysis
    dfd = df(idx.outlierStation); % drop the outliers
    sdif.thet.fconcs = NaN([size(dfd(1).concProf.gsClassModel.Cs), size(dfd,1)]);
    sdif.thet.MYconcs = NaN([size(dfd(1).concProf.MYfullPred.Cs), size(dfd,1)]);
    sdif.thet.cbis = NaN([1, size(dfd(1).concProf.gsClassModel.Cs, 2), size(dfd,1)]);
    for i = 1:size(dfd,1)
        sdif.thet.fconcs(:,:,i) = dfd(i).concProf.gsClassModel.Cs ./ 2650;
        sdif.thet.MYconcs(:,:,i) = dfd(i).concProf.MYfullPred.Cs ./ 2650;
        sdif.thet.cbis(1,:,i) = (dfd(i).concProf.gsClassModel.Cs(1,:) ./ 2650);
    end
    sdif.thet.Zs = linspace(0.05, 1, 51)';
    sdif.thet.msd = mean(sdif.thet.fconcs(3:end-2,:,:) - sdif.thet.MYconcs(3:end-2,:,:), 1); % mean signed deviation for every profile, when normalized in 0-1 cb space
    sdif.thet.nmsd = sdif.thet.msd ./ sdif.thet.cbis; % normalized to nearbed
    sdif.thet.sse = ( sum( (sdif.thet.fconcs - sdif.thet.MYconcs).^2, 1) / size(sdif.thet.MYconcs,1) ).^(0.5);
    sdif.thet.sdd = nanmean(sdif.thet.fconcs - sdif.thet.MYconcs, 3); %./  mean(cbis,3);
    sdif.thet.sddstd = std(sdif.thet.fconcs - sdif.thet.MYconcs, [], 3, 'omitnan');
    
    
    %% sediment diffusivity ratios
    stnsamples = arrayfun(@(x) size(x.waterSamplesTable,1), df, 'Unif', 1);
    nsamples = sum(stnsamples);
    sdif.tab.mat = NaN(nsamples*nClasses, 6);
    cnt = 1;
    settlingvelocitychart = get_DSV(df(1).waterSamplesTable(1,:).gsClass{:} .* 1e-6, 0.7, 3.5, con);
    for i = 1:size(df,1)
        for j = 1:stnsamples(i)
            for k = 1:nClasses
                sdif.tab.mat(cnt,1) = df(i).waterSamplesTable.sampleZ(j);
                sdif.tab.mat(cnt,2) = df(i).waterSamplesTable.sampleZnorm(j);
                sdif.tab.mat(cnt,3) = df(i).Velocity.ustarCalib / settlingvelocitychart(k);
                concNWbyClass = df(i).waterSamplesTable(j,:).concNWbyClass{1}(k);
                MYconcByClass = interp1(df(i).concProf.MYfullPred.Zs(:,k), df(i).concProf.MYfullPred.Cs(:,k), df(i).waterSamplesTable.sampleZ(j));
                sdif.tab.mat(cnt,4) = concNWbyClass;
                sdif.tab.mat(cnt,5) = MYconcByClass;
                sdif.tab.mat(cnt,6) = concNWbyClass / MYconcByClass;
                sdif.tab.mat(cnt,7) = df(i).waterSamplesTable(j,:).gsClass{1}(k);
                sdif.tab.mat(cnt,8) = year(df(i).CollectionDatetime);
                sdif.tab.mat(cnt,9) = (concNWbyClass-MYconcByClass)/concNWbyClass;
                cnt = cnt+1;
            end
        end
    end
    sdif.tab.tab = array2table(sdif.tab.mat,  'VariableNames', {'z', 'zh', 'ustarws', 'c', 'MY', 'ccMY', 'gs', 'cyear', 'pe'});
    
%     figure(); hold on;
%         set(gca,'ColorOrderIndex',1)
%         %% dashes for the fits
%         plot(df(1).concProf.gsClassModel.Cs, df(1).concProf.gsClassModel.Zs, '--')
%         plot(df(1).concProf.gsClassModel.CsSum, df(1).concProf.gsClassModel.ZsSum, 'k--')
%         set(gca,'ColorOrderIndex',1)
%         %% solids for the MY preds
%         plot(df(1).concProf.MYfullPred.Cs, df(1).concProf.gsClassPred.Zs, '-')
%         plot(df(1).concProf.MYfullPred.CsSum, df(1).concProf.gsClassPred.ZsSum, 'k-')
%     set(gca, 'xscale', 'log')
    
    
    
    
    
    
    
    
    
    
    
    %% make initial concentration examination plots
    if false
        f = figure(); hold on;
            plotYRdata(nearBed.ustarCalib ./ bed.d50ws, nearBed.concNW, idx.sampleYr, idx.sampleOut, f);
            xlabel('dimensionless shear velocity ($u_* / w_s$)')
            ylabel('near bed conc. (g/L)')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/dimlessUstarCalibVsnearBedconc.png');
                print('-depsc', '-painters', '-r300', './figsExport/dimlessUstarCalibVsnearBedconc.eps');
                pause(1);
            end

        f = figure(); hold on;
            mks1 = plotYRdata(nearBed.ustarCalib ./ bed.d10ws, nearBed.NWd10, idx.sampleYr, idx.sampleOut, f, 'Marker', 's');
            mks2 = plotYRdata(nearBed.ustarCalib ./ bed.d50ws, nearBed.NWd50, idx.sampleYr, idx.sampleOut, f, 'Marker', 'o');
            mk3 = plotYRdata(nearBed.ustarCalib ./ bed.d90ws, nearBed.NWd90, idx.sampleYr, idx.sampleOut, f, 'Marker', '^');
            set(gca, 'Xscale', 'log')
            xlabel('dimensionless shear velocity ($u_* / w_{s,\{D_{10},D_{50},D_{90}\}}$)')
            ylabel('grain size ($\mu m$)')
            legend([mks1(1), mks2(1), mk3(1)], {'$D_{10}$', '$D_{50}$', '$D_{90}$'}, 'Interpreter', 'latex')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/dimlessustarCalibVsgrainSize.png'); %#ok<*UNRCH>
                print('-depsc', '-painters', '-r300', './figsExport/dimlessUstarCalibVsgrainSize.eps');
                pause(1);
            end

        [ustar] = arrayfun(@(x) x.Velocity.ustarCalib, df);
        [~, sortustaridx] = sort(ustar);
        ustarinterpchart = linspace(min(ustar), max(ustar), 100);
        ustarcmap = parula(100); 
        f = figure(); hold on;
            for i = 1:size(df, 1)
                sortidx = sortustaridx(i);
                interpidx = find(ustar(sortidx) <= ustarinterpchart, 1, 'first');
                plot(gsDistClass2num(df(sortidx).bedData.gsDistBed), table2array(df(sortidx).bedData.gsDistBed), ...
                    'Color', ustarcmap(interpidx, :), 'LineWidth', 2)
                mat(:,i) = table2array(df(sortidx).bedData.gsDistBed);
            end
            plot(gsDistClass2num(df(sortidx).bedData.gsDistBed), mean(mat, 2), 'k-', 'LineWidth', 3)
            cb = colorbar('LimitsMode', 'manual', 'Limits', [min(ustar), max(ustar)]);
            cb.Label.String = '$u_*$ (m/s)';
            cb.Label.Interpreter = 'latex';
            caxis([min(ustar), max(ustar)])
            xlim([0 300])
            ylim([0 100])
            xlabel('grain size ($\mu m$)')
            ylabel('percent in class')
            title('channel bed gs dists, colored by $u_*$', 'interpreter', 'latex')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/grainSizeDistsVsColorUstar.png');
                print('-depsc', '-painters', '-r300', './figsExport/grainSizeDistsVsColorUstar.eps');
                pause(1);
            end

        f = figure(); hold on;
            for i = 1:size(df, 1)
                sortidx = sortustaridx(i);
                interpidx = find(ustar(sortidx) <= ustarinterpchart, 1, 'first');
                plot(gsDistClass2num(df(sortidx).bedData.gsDistBed), cumsum(table2array(df(sortidx).bedData.gsDistBed), 'omitnan'), ...
                    'Color', ustarcmap(interpidx, :), 'LineWidth', 2)
                mat(:,i) = table2array(df(sortidx).bedData.gsDistBed);
            end
            plot(gsDistClass2num(df(sortidx).bedData.gsDistBed), cumsum(mean(mat, 2), 'omitnan'), 'k-', 'LineWidth', 3)
            cb = colorbar('LimitsMode', 'manual', 'Limits', [min(ustar), max(ustar)]);
            cb.Label.String = '$u_*$ (m/s)';
            cb.Label.Interpreter = 'latex';
            caxis([min(ustar), max(ustar)])
            xlim([0 300])
            ylim([0 100])
            xlabel('grain size ($\mu m$)')
            ylabel('percent finer')
            title('channel bed gs dists, colored by u_{*}', 'interpreter', 'latex')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/grainSizeCumDistsVsColorUstar.png');
                print('-depsc', '-painters', '-r300', './figsExport/grainSizeCumDistsVsColorUstar.eps');
                pause(1);
            end
            
        f = figure(); hold on;
        depth_list = [0.95, 0.85, 0.75, 0.5, 0.1];
            for i = 1:length(depth_list)
                subplot(5, 1, i); hold on;
                [~, smpl_dims] = tableVarAtDepth(df, 'gsDistNW', depth_list(i));
                depthCells = cellfun(@(x) table2array(x), tableVarAtDepth(df, 'gsDistNW', depth_list(i)), 'unif', 0)';
                depthMat = [depthCells{:}];
                %sortidx = sortustaridx(i);
                %interpidx = find(ustar(sortidx) <= ustarinterpchart, 1, 'first');
                ustarvals = varRepLong(ustar, smpl_dims);
                text(0.8, 0.8, num2str(depth_list(i)), 'units', 'normalized');
                colormat = interp1(ustarinterpchart', ustarcmap, ustarvals);
                [h] = plot(gsDistClass2num(df(1).waterSamplesTable.gsDistNW{1}), depthMat, ...
                    'Color', 'k', 'LineWidth', 2);
                set(h, {'color'}, num2cell(colormat,2));
                set(gca, 'xscale', 'log', 'yscale', 'log')
                xlim([10, 300])
            end
            
            
         f = figure(); hold on;
            for i = 1:size(df, 1)
                sortidx = sortustaridx(i);
                interpidx = find(ustar(sortidx) <= ustarinterpchart, 1, 'first');
                ax1 = subplot(2, 1, 1); hold on;
                    plot(gsDistClass2num(df(sortidx).nearBedData.gsCharNearBedNWnorm), (table2array(df(sortidx).nearBedData.gsCharNearBedNWnorm)), ...
                        'Color', ustarcmap(interpidx, :), 'LineWidth', 2)
                    mat_char(:,i) = table2array(df(sortidx).nearBedData.gsCharNearBedNWnorm);
                ax2 = subplot(2, 1, 2); hold on;
                    plot(gsDistClass2num(df(sortidx).nearBedData.gsDistNearBedNWnorm), (table2array(df(sortidx).nearBedData.gsDistNearBedNWnorm)), ...
                        'Color', ustarcmap(interpidx, :), 'LineWidth', 2)
                    mat_dist(:,i) = table2array(df(sortidx).nearBedData.gsDistNearBedNWnorm);
            end
            subplot(2, 1, 1)
                plot(gsDistClass2num(df(sortidx).nearBedData.gsCharNearBedNWnorm), (mean(mat_char, 2)), 'k-', 'LineWidth', 3)
            subplot(2, 1, 2)
                plot(gsDistClass2num(df(sortidx).nearBedData.gsDistNearBedNWnorm), (mean(mat_dist, 2)), 'k-', 'LineWidth', 3)
            cb = colorbar(ax1, 'LimitsMode', 'manual', 'Limits', [min(ustar), max(ustar)]);
            cb.Label.String = '$u_*$ (m/s)';
            cb.Label.Interpreter = 'latex';
            cb = colorbar(ax2, 'LimitsMode', 'manual', 'Limits', [min(ustar), max(ustar)]);
            cb.Label.String = '$u_*$ (m/s)';
            cb.Label.Interpreter = 'latex';
            caxis(ax1, [min(ustar), max(ustar)])
            caxis(ax2, [min(ustar), max(ustar)])
            xlim(ax1, [0 300])
            ylim(ax1, [0 100])
            xlim(ax2, [0 300])
            ylim(ax2, [0 100])
            xlabel(ax2, 'grain size ($\mu m$)')
            ylabel(ax1, 'percent finer')
            ylabel(ax2, 'percent finer')
            title('near bed sample gs dists, colored by u_{*}', 'interpreter', 'latex')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/grainSizeCumDistsVsColorUstar.png');
                print('-depsc', '-painters', '-r300', './figsExport/grainSizeCumDistsVsColorUstar.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(nearBed.FlowDepth, nearBed.concWash, idx.sampleYr, idx.sampleOut, f);
            xlabel('flow depth (m)')
            ylabel('near bed conc. (g/L)')
            title('near-bed conc. w washload')
            box on
            axis square
            legend('2015', '2016', '2018')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/flowDepthVsnearBedconc.png');
                print('-depsc', '-painters', '-r300', './figsExport/flowDepthVsnearBedconc.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(nearBed.FlowDepth, nearBed.concNW, idx.sampleYr, idx.sampleOut, f);
            xlabel('flow depth (m)')
            ylabel('near bed conc. (g/L)')
            title('near-bed conc. w/o washload')
            box on
            axis square
            legend('2015', '2016', '2018')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/flowDepthVsnearBedconcNW.png');
                print('-depsc', '-painters', '-r300', './figsExport/flowDepthVsnearBedconcNW.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(nearBed.FlowDepth, nearBed.ustarCalib, idx.sampleYr, idx.sampleOut, f);
            xlabel('flow depth (m)')
            ylabel('$u_*$ (m/s)')
            title('')
            box on
            axis square
            legend('2015', '2016', '2018')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/flowDepthVsUstar.conc.png');
                print('-depsc', '-painters', '-r300', './figsExport/flowDepthVsUstar.conc.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(nearBed.ustarCalib ./ bed.d50ws .* nearBed.Rep .^ 0.6, nearBed.concNW, idx.sampleYr, idx.sampleOut, f);
            xlabel('$(u_* / w_s) \textrm{Rep}^{0.6}$')
            ylabel('near bed conc. (g/L)')
            title('')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/dimlessUstarRepVsnearBed.concNW.png');
                print('-depsc', '-painters', '-r300', './figsExport/dimlessUstarRepVsnearBed.concNW.eps');
                pause(1);
            end
    end
        
    
    %% 2016/2018 DSP and Lotw fit matches
    % attempt fit for discharge and depth to explain ustar better
    velocityCalibTable.logDischarge = log10(velocityCalibTable.discharge);
    velocityCalibTable.logDepth = log10(velocityCalibTable.depth);
    velocityCalibTable.logUstarFric = log10(velocityCalibTable.ustarFric);
    velocityCalibTable.depthSq= velocityCalibTable.depth .* velocityCalibTable.depth;
    velocityCalibTable.dischargeOverDepthSq= velocityCalibTable.discharge ./ (velocityCalibTable.depth .* velocityCalibTable.depth);
    velocityCalibTable.dischargeNormed= velocityCalibTable.discharge ./ 3000;
    velocityCalibTable.ustarNoAlphaSquared = velocityCalibTable.ustarNoAlpha .* velocityCalibTable.ustarNoAlpha;
    velocityCalibTable.ustarFricSquared = velocityCalibTable.ustarFric .* velocityCalibTable.ustarFric;
    velocityCalibTable.ustarFricTilde = velocityCalibTable.ustarFric / ((con.R * con.g * con.nu)^(1/3));
    velocityCalibTable.depthTilde = (velocityCalibTable.depth * (con.g^(1/3))) / ((con.R*con.nu)^(2/3));
    velocityCalibTable.dischargeTilde = velocityCalibTable.discharge ./ (sqrt(con.g .* ((90e-6)^5)));
    
    lmAlpha = fitlm(velocityCalibTable,'ustarAlpha ~ logDischarge + depth');
    lmNoAlpha = fitlm(velocityCalibTable,'ustarNoAlpha ~ logDischarge + depth');
    lmDSP = fitlm(velocityCalibTable,'ustarDSP ~ logDischarge + depth');
    lmFric = fitlm(velocityCalibTable,'ustarFric ~ logDischarge + depth');
    lmFricHomo = fitlm(velocityCalibTable,'ustarFric ~ dischargeOverDepthSq');
    lmFricHomo2 = fitlm(velocityCalibTable,'ustarFric ~ dischargeOverDepthSq + ustarDSP');
    lmLoggedDimless = fitlm(velocityCalibTable,'logUstarFric ~ logDischarge + logDepth');
    lmAlphafunc = @(Q,H) lmAlpha.Coefficients.Estimate(3).*log10(Q) + lmAlpha.Coefficients.Estimate(2).*H + lmAlpha.Coefficients.Estimate(1);
    lmNoAlphafunc = @(Q,H) lmNoAlpha.Coefficients.Estimate(3).*log10(Q) + lmNoAlpha.Coefficients.Estimate(2).*H + lmNoAlpha.Coefficients.Estimate(1);
    lmDSPfunc = @(Q,H) lmDSP.Coefficients.Estimate(3).*log10(Q) + lmDSP.Coefficients.Estimate(2).*H + lmDSP.Coefficients.Estimate(1);
    lmFricfunc = @(Q,H) lmFric.Coefficients.Estimate(3).*log10(Q) + lmFric.Coefficients.Estimate(2).*H + lmFric.Coefficients.Estimate(1);
    lmFricHomofunc = @(Q,H) lmFricHomo.Coefficients.Estimate(2) .* (Q./(H.*H)) + lmFricHomo.Coefficients.Estimate(1);
    lmFricHomo2func = @(Q,H) lmFricHomo2.Coefficients.Estimate(2) .* (sqrt(9.81*H*6.4e-5)) + lmFricHomo2.Coefficients.Estimate(3) .* (Q./(H.*H)) ...
        + lmFricHomo2.Coefficients.Estimate(1);
    lmLoggedDimlessfunc = @(Q,H) lmLoggedDimless.Coefficients.Estimate(2).*log10(Q) + lmLoggedDimless.Coefficients.Estimate(3).*log10(H) + lmLoggedDimless.Coefficients.Estimate(1);
    eAlpha = lmAlphafunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    eNoAlpha = lmNoAlphafunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    eDSP = lmDSPfunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    eFric = lmFricfunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    eFricHomo = lmFricHomofunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    eFricHomo2 = lmFricHomo2func(velocityCalibTable.discharge, velocityCalibTable.depth);
    eLoggedDimless = lmLoggedDimlessfunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    % nonlinear model below here
    nlmNoAlpha = fitnlm(velocityCalibTable,'ustarNoAlpha ~ a*(logDischarge^m) * (depth^n)',[1 1 1]);
    nlmNoAlphafunc = @(Q,H) nlmNoAlpha.Coefficients.Estimate(1).*(log10(Q)).^nlmNoAlpha.Coefficients.Estimate(2) .* ...
            H.^nlmNoAlpha.Coefficients.Estimate(3);
    enlmNoAlpha = nlmNoAlphafunc(velocityCalibTable.logDischarge, velocityCalibTable.depth);
    nlmFric = fitnlm(velocityCalibTable,'ustarFric ~ a*(logDischarge^m) * (depth^n)',[0 1 1]);
    nlmFricfunc = @(Q,H) nlmFric.Coefficients.Estimate(1).*(log10(Q)).^nlmFric.Coefficients.Estimate(2) .* ...
            H.^nlmFric.Coefficients.Estimate(3);
    enlmFric = nlmFricfunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    nlmFricHomo = fitnlm(velocityCalibTable,'ustarFricSquared ~ ((a* discharge ^ b) / (400 * depth)) * (9.81*0.000064)',[1 1]);
    nlmFricHomofunc = @(Q,H) ((nlmFricHomo.Coefficients.Estimate(1) .* Q .^ nlmFricHomo.Coefficients.Estimate(2)) ./ 400 ./ (H)) .* (9.81*0.000064);
    enlmFricHomo = nlmFricHomofunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    nlmFricHomoH2 = fitnlm(velocityCalibTable,'ustarFric ~ a*(discharge / (depth^2))^n',[0 1]);
    nlmFricHomoH2func = @(Q,H) nlmFricHomoH2.Coefficients.Estimate(1) .* (Q ./ H.^2) .^ nlmFricHomoH2.Coefficients.Estimate(2);
    enlmFricHomoH2 = nlmFricHomoH2func(velocityCalibTable.discharge, velocityCalibTable.depth);
    nlmFricHomoDSP = fitnlm(velocityCalibTable,'ustarFricSquared ~ (a * (discharge/3000) ^ b) * (depth*9.81*0.000064)',[1 1]);
    nlmFricHomoDSPfunc = @(Q,H) (nlmFricHomoDSP.Coefficients.Estimate(1) .* (Q / 3000) .^ nlmFricHomoDSP.Coefficients.Estimate(2)) .* (H.*9.81*0.000064);
    enlmFricHomoDSP = nlmFricHomoDSPfunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    lmFricHomoDSP = fitlm(velocityCalibTable,'ustarFric ~ dischargeNormed + ustarDSP');
    lmFricHomoDSPfunc = @(Q,H) (lmFricHomoDSP.Coefficients.Estimate(3) .* (Q / 3000) + lmFricHomoDSP.Coefficients.Estimate(2) .* sqrt(H.*9.81*0.000064) + lmFricHomoDSP.Coefficients.Estimate(1));
    elmFricHomoDSP = lmFricHomoDSPfunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    nlmFricHomoTilde = fitnlm(velocityCalibTable,'ustarFricTilde ~ a*(dischargeTilde * depthTilde)^n', [1 1]);
    nlmFricHomoTildefunc = @(Q,H) nlmFricHomoTilde.Coefficients.Estimate(1) .* ((Q./sqrt(9.81.*(90e-6).^5)) .* ((H.*(con.g^(1/3)))/((con.R*con.nu)^(2/3)))) .^ nlmFricHomoTilde.Coefficients.Estimate(2);
    enlmFricHomoTilde = nlmFricHomoTildefunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    lmFricHomoTilde = fitlm(velocityCalibTable,'ustarFricTilde ~ dischargeTilde + depthTilde');
    lmFricHomoTildefunc = @(Q,H) lmFricHomoTilde.Coefficients.Estimate(3) .* ((Q./sqrt(con.g*(90e-6)^5)) + lmFricHomoTilde.Coefficients.Estimate(2) .* ((H.*(con.g^(1/3)))/((con.R*con.nu)^(2/3)))) + lmFricHomoTilde.Coefficients.Estimate(1);
    elmFricHomoTilde = lmFricHomoTildefunc(velocityCalibTable.discharge, velocityCalibTable.depth);
    ustarUnTildeFunc = @(ustarTilde) ustarTilde .* (con.R .* con.g .* con.nu).^(1/3);
    % single coefficient model to DSP
    lmSingle = fitnlm(velocityCalibTable,'ustarFric ~ (9.81 * depth * 0.000064 * A)', [1]);
    lmSinglefunc = @(H, S) sqrt(con.g .* H .* lmSingle.Coefficients.Estimate(1) .* S);
    eSingle = lmSinglefunc(velocityCalibTable.depth, velocityCalibTable.slope);
        
    % choose to save one of the models as the calibration and update the data in the table
%     ustar_calib_model = nlmFric;
%     ustar_calib_func = nlmFricfunc;
%     ustar_calib_func_str = [num2str(round(nlmFric.Coefficients.Estimate(1), 4)), 'log10(Q) ^{', num2str(round(nlmFric.Coefficients.Estimate(2),3)), '}H ^{', num2str(round(nlmFric.Coefficients.Estimate(3),3)) '}'];
    ustar_calib_model = lmLoggedDimless;
    ustar_calib_func = @(Q,H) 10^lmLoggedDimless.Coefficients.Estimate(1) .* (Q.^lmLoggedDimless.Coefficients.Estimate(2)) .* (H.^lmLoggedDimless.Coefficients.Estimate(3));
    ustar_calib_func_str = [num2str(round(10^lmLoggedDimless.Coefficients.Estimate(1), 4)), '(Q^{', num2str(round(lmLoggedDimless.Coefficients.Estimate(2),3)), '}H ^{', num2str(round(lmLoggedDimless.Coefficients.Estimate(3),3)) '}'];
    save('./dataExport/ustar_calibration.mat', 'ustar_calib_model', 'ustar_calib_func', 'ustar_calib_func_str')
    velocityCalibTable.ustarCalib = ustar_calib_func(velocityCalibTable.discharge, velocityCalibTable.depth);
    station.velCalib.ustarIntoCalib = station.velCalib.ustarFric;
    station.velCalib.ustarOutofCalib = velocityCalibTable.ustarCalib;
    
    % COMPARE REGRESSIONS
    eeLoggedDimless = 10.^(eLoggedDimless);
    Bbar = mean(enlmFric);
    SStot = sum((enlmFric - Bbar).^2);
    SSreg = sum((eeLoggedDimless - Bbar).^2);
    SSres = sum((enlmFric - eeLoggedDimless).^2);
    R2 = 1 - SSres/SStot;
    figure()
    subplot(3, 3, 1); hold on;
        plot(velocityCalibTable.ustarFric, enlmFric, 'o')
        title(['$a (\log_{10} Q_w) ^{b} H ^{c}$'], 'FontSize', 14)
        plot([0, 0.1], [0, 0.1], 'k--')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex'); box on;
        text(0.7, 0.1, ['$R^2 = $', num2str(nlmFric.Rsquared.Ordinary, 2)], 'Units', 'normalized', 'interpreter', 'latex')
        text(0.1, 0.9, 'a', 'Units', 'normalized', 'interpreter', 'latex')
    subplot(3, 3, 2); hold on;
        plot(velocityCalibTable.ustarFric, enlmFricHomoH2, 'o')
        title('${u_*} = a (Q_w / H^2) ^b$', 'FontSize', 14)
        plot([0, 0.1], [0, 0.1], 'k--')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex'); box on;
        text(0.7, 0.1, ['$R^2 = $', num2str(nlmFricHomoH2.Rsquared.Ordinary, 2)], 'Units', 'normalized', 'interpreter', 'latex')
        text(0.1, 0.9, 'b', 'Units', 'normalized', 'interpreter', 'latex')
    subplot(3, 3, 3); hold on;
        plot(velocityCalibTable.ustarFric, sqrt(enlmFricHomoDSP), 'o')
        title('${u_*}^2 = a (Q_w/Q_{w,bf}) ^ b (gHS)$', 'FontSize', 14)
        plot([0, 0.1], [0, 0.1], 'k--')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex'); box on;
        text(0.7, 0.1, ['$R^2 = $', num2str(nlmFricHomoDSP.Rsquared.Ordinary, 2)], 'Units', 'normalized', 'interpreter', 'latex')
        text(0.1, 0.9, 'c', 'Units', 'normalized', 'interpreter', 'latex')        
    subplot(3, 3, 4); hold on;
        plot(velocityCalibTable.ustarFric, sqrt(enlmFricHomo), 'o')
        title('${u_*}^2 = \frac{a{Q_w}^b}{B H} gS$', 'FontSize', 14)
        plot([0, 0.1], [0, 0.1], 'k--')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex'); box on;
        text(0.7, 0.1, ['$R^2 = $', num2str(nlmFricHomo.Rsquared.Ordinary, 2)], 'Units', 'normalized', 'interpreter', 'latex')
        text(0.1, 0.9, 'd', 'Units', 'normalized', 'interpreter', 'latex')
    subplot(3, 3, 5); hold on;
        plot(velocityCalibTable.ustarFric, elmFricHomoDSP, 'o')
        title('${u_*} = a(Q_w/Q_{w,bf}) + b\sqrt{gHS} + c$', 'FontSize', 14)
        plot([0, 0.1], [0, 0.1], 'k--')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex'); box on;
        text(0.7, 0.1, ['$R^2 = $', num2str(lmFricHomoDSP.Rsquared.Ordinary, 2)], 'Units', 'normalized', 'interpreter', 'latex')
        text(0.1, 0.9, 'e', 'Units', 'normalized', 'interpreter', 'latex')
    subplot(3, 3, 6); hold on;
        plot(velocityCalibTable.ustarFric, ustarUnTildeFunc(enlmFricHomoTilde), 'o')
        title('${u_*} =  a ( \tilde{Q_w} \tilde{H} ) ^b$', 'FontSize', 14)
        plot([0, 0.1], [0, 0.1], 'k--')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex'); box on;
        text(0.7, 0.1, ['$R^2 = $', num2str(nlmFricHomoTilde.Rsquared.Ordinary, 2)], 'Units', 'normalized', 'interpreter', 'latex')
        text(0.1, 0.9, 'f', 'Units', 'normalized', 'interpreter', 'latex')
    subplot(3, 3, 7); hold on;
        plot(velocityCalibTable.ustarFric, ustarUnTildeFunc(elmFricHomoTilde), 'o')
        title('$\tilde{u_*} = a \tilde{Q_w} + b \tilde{H} + c$', 'FontSize', 14)
        plot([0, 0.1], [0, 0.1], 'k--')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex'); box on;
        text(0.7, 0.1, ['$R^2 = $', num2str(lmFricHomoTilde.Rsquared.Ordinary, 2)], 'Units', 'normalized', 'interpreter', 'latex')
        text(0.1, 0.9, 'g', 'Units', 'normalized', 'interpreter', 'latex')
    subplot(3, 3, 8); hold on;
        plot(velocityCalibTable.ustarFric, 10.^(eLoggedDimless), 'o')
        title('$\log_{10}{u_*} = b_0 + b_1 \log_{10}{Q_w} + b_2 \log_{10}H$', 'FontSize', 14)
        plot([0, 0.1], [0, 0.1], 'k--')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex'); box on;
        text(0.7, 0.1, ['$R^2 = $', num2str(lmLoggedDimless.Rsquared.Ordinary, 2)], 'Units', 'normalized', 'interpreter', 'latex')
        text(0.1, 0.9, 'h', 'Units', 'normalized', 'interpreter', 'latex')
    subplot(3, 3, 9); hold on;
        plot(enlmFric, eeLoggedDimless, 'o')
        plot([0, 0.1], [0, 0.1], 'k--')
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex'); box on;
        text(0.7, 0.1, ['$R^2 = $', num2str(R2, 4)], 'Units', 'normalized', 'interpreter', 'latex')
        text(0.1, 0.9, 'i', 'Units', 'normalized', 'interpreter', 'latex')
        xlabel('old values (a)')
        ylabel('new values (h)')
        
        
    j=1;
    if false
        f = figure(); hold on;
            plotYRdata(station.ustarDSP, station.ustarCalib, idx.stationYr, idx.stationOut, f)
            plot([0 0.15], [0 0.15], 'k-', 'LineWidth', 1);
            axis square
            xlabel('$u_*$ from true DSP (m/s)')
            ylabel('$u_*$ from calibration (m/s)')
            box on
    %         legend('2016', '2018', 'Location', 'NorthWest')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/DSPandFitCorrelation.png');
                print('-depsc', '-painters', '-r300', './figsExport/DSPandFitCorrelation.eps');
                pause(1);
            end

        f = figure(); hold on;
            plot(station.velCalib.ustarNoAlpha(station.velCalib.ustarYr==2016), eNoAlpha(station.velCalib.ustarYr==2016), ...
                'Marker', mar2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col2, 'MarkerSize', 6, 'LineStyle', 'none')
            plot(station.velCalib.ustarNoAlpha(station.velCalib.ustarYr==2018), eNoAlpha(station.velCalib.ustarYr==2018), ... 
                'Marker', mar3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col3, 'MarkerSize', 6, 'LineStyle', 'none')
            plot([0 0.15], [0 0.15], 'k-', 'LineWidth', 1);
            axis square
            xlabel('$u_*$ from profile data fit NoAlpha (m/s)')
            ylabel('$u_*$ from coeff on DSP calibration (m/s)')
            box on
            legend('2016', '2018', 'Location', 'NorthWest')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/ustarNoAlphaAndCalibrationCorrelation.png');
                print('-depsc', '-painters', '-r300', './figsExport/ustarNoAlphaAndCalibrationCorrelation.eps');
                pause(1);
            end

        f = figure(); hold on;
            plot(station.velCalib.ustarAlpha(station.velCalib.ustarYr==2016), eAlpha(station.velCalib.ustarYr==2016), ...
                'Marker', mar2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col2, 'MarkerSize', 6, 'LineStyle', 'none')
            plot(station.velCalib.ustarAlpha(station.velCalib.ustarYr==2018), eAlpha(station.velCalib.ustarYr==2018), ... 
                'Marker', mar3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col3, 'MarkerSize', 6, 'LineStyle', 'none')
            plot([0 0.15], [0 0.15], 'k-', 'LineWidth', 1);
            axis square
            xlabel('$u_*$ from profile data fit with Alpha (m/s)')
            ylabel('$u_*$ from coeff on DSP calibration (m/s)')
            box on
            legend('2016', '2018', 'Location', 'NorthWest')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/ustarAlphaAndCalibrationCorrelation.png');
                print('-depsc', '-painters', '-r300', './figsExport/ustarAlphaAndCalibrationCorrelation.eps');
                pause(1);
            end

        f = figure(); hold on;
            plot(station.velCalib.ustarDSP(station.velCalib.ustarYr==2016), eDSP(station.velCalib.ustarYr==2016), ...
                'Marker', mar2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col2, 'MarkerSize', 6, 'LineStyle', 'none')
            plot(station.velCalib.ustarDSP(station.velCalib.ustarYr==2018), eDSP(station.velCalib.ustarYr==2018), ... 
                'Marker', mar3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col3, 'MarkerSize', 6, 'LineStyle', 'none')
            plot([0 0.15], [0 0.15], 'k-', 'LineWidth', 1);
            axis square
            xlabel('$u_*$ from true DSP data (m/s)')
            ylabel('$u_*$ from coeff on DSP calibration (m/s)')
            box on
            legend('2016', '2018', 'Location', 'NorthWest')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/ustarDSPAndCalibrationCorrelation.png');
                print('-depsc', '-painters', '-r300', './figsExport/ustarDSPAndCalibrationCorrelation.eps');
                pause(1);
            end

        f = figure(); hold on;
            plot(station.velCalib.discharge(station.velCalib.ustarYr==2016), station.velCalib.depth(station.velCalib.ustarYr==2016), ...
                'Marker', mar2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col2, 'MarkerSize', 6, 'LineStyle', 'none')
            plot(station.velCalib.discharge(station.velCalib.ustarYr==2018), station.velCalib.depth(station.velCalib.ustarYr==2018), ...
                'Marker', mar3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col3, 'MarkerSize', 6, 'LineStyle', 'none')
            axis square
            xlabel('discharge at Lijin (m$^3$/s)')
            ylabel('depth at measurement (m)')
            box on
            legend('2016', '2018', 'Location', 'SouthEast')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/dischargeVsDepthCorrelation.png');
                print('-depsc', '-painters', '-r300', './figsExport/dischargeVsDepthCorrelation.eps');
                pause(1);
            end

        f = figure(); hold on;
            plot(station.velCalib.discharge(station.velCalib.ustarYr==2016), station.velCalib.ustarIntoCalib(station.velCalib.ustarYr==2016), ...
                'Marker', mar2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col2, 'MarkerSize', 6, 'LineStyle', 'none')
            plot(station.velCalib.discharge(station.velCalib.ustarYr==2018), station.velCalib.ustarIntoCalib(station.velCalib.ustarYr==2018), ...
                'Marker', mar3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col3, 'MarkerSize', 6, 'LineStyle', 'none')
            axis square
            xlabel('discharge at Lijin (m$^3$/s)')
            ylabel('$u_*$ that went into calib (m/s)')
            box on
            legend('2016', '2018', 'Location', 'NorthWest')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/dischargeVsUstarIntoCalibCorrelation.png');
                print('-depsc', '-painters', '-r300', './figsExport/dischargeVsUstarIntoCalibCorrelation.eps');
                pause(1);
            end

        f = figure(); hold on;
            plot(station.velCalib.depth(station.velCalib.ustarYr==2016), station.velCalib.ustarIntoCalib(station.velCalib.ustarYr==2016), ...
                'Marker', mar2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col2, 'MarkerSize', 6, 'LineStyle', 'none')
            plot(station.velCalib.depth(station.velCalib.ustarYr==2018), station.velCalib.ustarIntoCalib(station.velCalib.ustarYr==2018), ...
                'Marker', mar3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col3, 'MarkerSize', 6, 'LineStyle', 'none')
            [stdDSP] = plot(0:0.2:8, (sqrt(con.g .* (0:0.2:8) .* 6.4e-5)), 'k--');
            [stdDSP] = plot(0:0.2:8, (sqrt(con.g .* (0:0.2:8) .* 1e-4)), 'k--');
            [stdDSP] = plot(0:0.2:8, (sqrt(con.g .* (0:0.2:8) .* 3e-4)), 'k--');
            text(5.5, 0.12, '$3 \times 10^{-4}$', 'rotation', 21)
            text(5.5, 0.085, '$1 \times 10^{-4}$', 'rotation', 14)
            text(5.5, 0.0685, '$6 \times 10^{-5}$', 'rotation', 11.5)
            axis square
            xlabel('flow depth (m)')
            ylabel('$u_*$ measured from profiles (m/s)')
            box on
            legend('2016', '2018', 'DSP', 'Location', 'SouthEast')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/depthVsUstarIntoCalibCorrelation.png');
                print('-depsc', '-painters', '-r300', './figsExport/depthVsUstarIntoCalibCorrelation.eps');
                pause(1);
            end


        f = figure(); hold on; % this plot should be the one to demonstrate that the calibration shifts away from the DSP: includes all data and both trend lines
            [mkData] = plotYRdata(station.depth, station.ustarCalib, idx.stationYr, idx.stationOut, f);
            [stdDSP] = plot(0:0.2:8, (sqrt(con.g .* (0:0.2:8) .* df(1).Velocity.slope)), 'k--');
    %         axis square
            xlabel('depth (m)')
            ylabel('$u_*$ (m/s)')
            box on
            legend([mkData(1), stdDSP], {'data after calibration', 'std DSP'}, 'Location', 'SouthEast')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/depthVsUstarCalibration.png');
                print('-depsc', '-painters', '-r300', './figsExport/depthVsUstarCalibration.eps');
                pause(1);
            end
    end
        
    %% all velocity profiles
    velProfIdxLocs = find(velProfIdx);
    if false
        figure();
        for p = 1:sum(velProfIdx)
            subplot(3, 5, p); hold on;
            dfIdx = velProfIdxLocs(p);
            if datenum(df(dfIdx).CollectionDate) < datenum('01/01/2016') 
                yrCol = col1;
                yrMark = mar1;
            elseif datenum(df(dfIdx).CollectionDate) < datenum('01/01/2017')
                yrCol = col2;
                yrMark = mar2;
            else
                yrCol = col3;
                yrMark = mar3;
            end
            if ~isempty(df(dfIdx).Velocity.adcp.mean)
                errorbar(df(dfIdx).Velocity.adcpTable.mean, df(dfIdx).Velocity.adcpTable.measZ, df(dfIdx).Velocity.adcpTable.std, 'k.', 'horizontal')
                plot(df(dfIdx).Velocity.adcpTable.mean, df(dfIdx).Velocity.adcpTable.measZ, 'k', 'MarkerFaceColor', yrCol, 'MarkerSize', 5, 'Marker', yrMark, 'LineStyle', 'none')
                plot(df(dfIdx).velProf.modelNoAlpha.Us, df(dfIdx).velProf.modelNoAlpha.Zs, 'k-')
            end
            if ~isempty(df(dfIdx).Velocity.meter.velmean)
                if ~iscell(df(dfIdx).Velocity.meterTable.err(1))
                    errorbar(df(dfIdx).Velocity.meterTable.mean, df(dfIdx).Velocity.meterTable.measZ, df(dfIdx).Velocity.meterTable.err, 'k.', 'horizontal')
                else
                    df(dfIdx).Velocity.meterTable.errneg = df(dfIdx).Velocity.meterTable.mean - arrayfun(@(x) x{1}(1), df(dfIdx).Velocity.meterTable.err, 'Unif', 1);
                    df(dfIdx).Velocity.meterTable.errpos = arrayfun(@(x) x{1}(2), df(dfIdx).Velocity.meterTable.err, 'Unif', 1) - df(dfIdx).Velocity.meterTable.mean;
                    errorbar(df(dfIdx).Velocity.meterTable.mean, df(dfIdx).Velocity.meterTable.measZ, df(dfIdx).Velocity.meterTable.errneg, df(dfIdx).Velocity.meterTable.errpos, 'k.', 'horizontal')
                end
                plot(df(dfIdx).Velocity.meterTable.mean, df(dfIdx).Velocity.meterTable.measZ, 'k', 'MarkerSize', 5, 'Marker', yrMark, 'LineStyle', 'none')
                plot(df(dfIdx).velProf.modelNoAlpha.Us, df(dfIdx).velProf.modelNoAlpha.Zs, 'k-')
            end
            WPvel_mdl = plot(df(dfIdx).velProf.WPfullPred.Us, df(dfIdx).concProf.WPfullPred.Zs, ...
                    'r--', 'LineWidth', 1.5); % wright parker prediction
            params1a = {strcat('$z_0 = ', sprintsci(df(dfIdx).velProf.modelNoAlpha.params.z0), '$'), ...
            strcat('$u_* = ', num2str(round(df(dfIdx).velProf.modelNoAlpha.params.ustar,4)), '$'), ...
            strcat('$R^2 = ', num2str(round(df(dfIdx).velProf.modelNoAlpha.params.model.Rsquared.Ordinary,2)), '$')};
            text(0.05, 0.7, sprintf('%s\n', params1a{:}), 'Parent', gca, 'units', 'normalized');
            plot(xlim, repmat(df(dfIdx).FlowDepthAtCollection, 1, 2), 'Color', [0 0 1], 'LineStyle', '--', 'LineWidth', 1.5);
            ylim([0 df(dfIdx).FlowDepthAtCollection*1.1])
            xlabel('velocity (m/s)')
            ylabel('distance above bed (m)')
            title(df(dfIdx).StationID, 'interpreter', 'none')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 1400 1000], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
        end
        if printOut
            print('-dpng', '-r300', './figsExport/everyVelProf.png');
            print('-depsc', '-painters', '-r300', './figsExport/everyVelProf.eps');
            pause(1);
        end
    end
          
    %% flow depth grouped concentrations vs shear velocity regressions
    regConcNWCell = arrayfun(@(x) x.waterSamplesTable.concNW, df, 'Unif', 0);
    regConcNW = vertcat(regConcNWCell{:});
    regZnormCell = arrayfun(@(x) x.waterSamplesTable.sampleZnorm, df, 'Unif', 0);
    regZnorm = vertcat(regZnormCell{:});
    regVarDims = cellfun(@length, regConcNWCell);
    [regUstarDSP] = varVsDepth(df, 'Velocity.ustarDSP', regVarDims);
    % figure()
    % gscatter(regUstarDSP, regConcNW, categorical(regZnorm))
    
    %% station KT1-C through the 2015 flood
    stationLoc.location = 'KT1-C';
    stationLoc.Idx = logical(  varRepLong( double(arrayfun(@(x) strcmp(x.StationLocation, stationLoc.location), df)), nearBed.dims )  );
    stationLoc.UstarDSP = nearBed.ustarCalib(stationLoc.Idx);
    stationLoc.NearBed.NWd50ws = bed.d50ws(stationLoc.Idx);
    stationLoc.nearBed.concNW = nearBed.concNW(stationLoc.Idx);
    if false
        figure()
            plot(stationLoc.UstarDSP./stationLoc.NearBed.NWd50ws, stationLoc.nearBed.concNW, 'ks', 'MarkerFaceColor', col1, 'MarkerSize', 6)
            xlabel('dimensionless shear velocity $(u_* / w_s)$')
            ylabel('near bed conc. (g/L)')
            title('KT1-C only')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
    end

    %% transect KT1-X through the 2015 flood
    transectLoc.transect = 'KT1';
    transectLoc.idx = logical(  varRepLong( double(arrayfun(@(x) strcmp(x.Transect, transectLoc.transect), df)), nearBed.dims )  );
    transectLoc.ustarCalib = nearBed.ustarCalib(transectLoc.idx);
    transectLoc.nearBed.NWd50ws = bed.d50ws(transectLoc.idx);
    transectLoc.nearBed.concNW = nearBed.concNW(transectLoc.idx);
    if false
        figure()
            plot(transectLoc.ustarCalib./transectLoc.nearBed.NWd50ws, transectLoc.nearBed.concNW, 'ks', 'MarkerFaceColor', col1, 'MarkerSize', 6)
            xlabel('dimensionless shear velocity $(u_* / w_s)$')
            ylabel('near bed conc. (g/L)')
            title('KT1 only')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
    end
    
    %% concProf plots -- theta
    if false
        f = figure(); hold on;  % plot the theta norm as a function of depth
            for i = 1:1:5
                subplot(5, 1, (5-i+1))
                mean_range = (i*10-9):(i*10-1);
                plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted.*station.suspension.bedDistWeightedRep.^(0.6), ...
                    nanmean(station.theta.byDepthToGsClassModelNorm(:, mean_range), 2), idx.stationYr, idx.stationOut, f);
                plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
                ylabel('$\Theta / c_b$')
                text(0.1, 0.1, ['depth=',num2str((5-i+1)/5)], 'Units', 'Normalized')
                ylim([0, 10])
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            end
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsThetaNormByDepthToTotalModelAsDepth.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsThetaNormByDepthToTotalModelAsDepth.eps');
                pause(1);
            end
        
        f = figure(); hold on;  % plot the theta norm as a function of depth for WP models
            for i = 1:1:5
                subplot(5, 1, (5-i+1))
                mean_range = (i*10-9):(i*10-1);
                plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted.*station.suspension.bedDistWeightedRep.^(0.6), ...
                    nanmean(station.theta.byDepthWPfullGsClassToGsClassModelNorm(:, mean_range), 2), ...
                    idx.stationYr, idx.stationOut, f);
                plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
                ylabel('$\Theta / c_b$')
                text(0.1, 0.1, ['depth=',num2str((5-i+1)/5)], 'Units', 'Normalized')
                ylim([-5, 5])
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            end
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsThetaNormbyDepthWPfullGsClassToGsClassModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsThetaNormbyDepthWPfullGsClassToGsClassModel.eps');
                pause(1);
            end
        
        f = figure(); hold on;  % plot the theta norm as a function of depth and grain size for WP models
        m_plot = 0.015; % margins
        s_plot = 0.015; % spacing
        nx = nClasses;
        ny = 5;
        w_plot = (1-2*m_plot-(nx-1)*s_plot) / nx;
        h_plot = (1-2*m_plot-(ny-1)*s_plot) / ny;
        x_pos = repmat(m_plot+(0:(nx-1))*w_plot+(0:(nx-1))*s_plot, ny, 1)';
        y_pos = repmat(1-(h_plot+m_plot+(0:(ny-1))*h_plot+(0:(ny-1))*s_plot), nx, 1);
        [gs_mesh, z_mesh] = meshgrid(1:nClasses, 1:5);
            for i = 1:5
                mean_range = (i*10-9):(i*10-1);
                for j = 1:nClasses
                    % subplot(5, nClasses, i*(j-1)+j)
                    subplot('Position', [x_pos(j, i), y_pos(j, i), w_plot , h_plot])
                    
                    plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted.*station.suspension.bedDistWeightedRep.^(0.6), ...
                        squeeze( nanmean(  station.theta.byDepthWPfullGsClassToGsClassModelDepthNorm(mean_range, j, :)  , 1) ), ...
                        idx.stationYr, idx.stationOut, f);
                    plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
                    % ylabel('$\Theta / c_b$')
                    text(0.1, 0.1, ['depth=',num2str((5-i+1)/5)], 'Units', 'Normalized')
                    text(0.1, 0.2, ['gs=',num2str(j)], 'Units', 'Normalized')
                    ylim([-1, 1])
                    box on
                    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
                end
            end
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsThetaNormbyDepthWPfullGsClassToGsClassModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsThetaNormbyDepthWPfullGsClassToGsClassModel.eps');
                pause(1);
            end
            

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormD50, station.theta.meanNormToTotalModel, idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{50},nb})$')
            ylabel('$\Theta_{t} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsThetaNormToTotalModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsThetaNormToTotalModel.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted, station.theta.meanNormToTotalModel, idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},nb})$')
            ylabel('$\Theta_{t} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormDistWeightedVsThetaNormToTotalModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormDistWeightedVsThetaNormToTotalModel.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberBedD50, station.theta.meanNormToTotalModel, idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{50},b})$')
            ylabel('$\Theta_{t} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberBedD50VsThetaNormToTotalModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberBedD50VsThetaNormToTotalModel.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberBedDistWeighted, station.theta.meanNormToTotalModel, idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},b})$')
            ylabel('$\Theta_{t} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberBedDistWeightedVsThetaNormToTotalModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberBedDistWeightedVsThetaNormToTotalModel.eps');
                pause(1);
            end

         f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted .* station.suspension.nearBedNWnormDistWeightedRep .^ 0.6, station.theta.meanNormToTotalModel, ...
                idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $\times$ Rep $([u_* / w_{s,D_{w},nb}] \cdot $Rep$^{0.6})$')
            ylabel('$\Theta_{t} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormRepD50VsThetaNormToTotalModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormRepD50VsThetaNormToTotalModel.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.ustarCalib, station.theta.meanNormToTotalModel, idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('shear velocity $(u_*)$')
            ylabel('$\Theta_{t} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/UstarVsThetaNormToTotalModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/UstarVsThetaNormToTotalModel.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormD50, station.theta.meanNormToGsClassModel, idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{50},nb})$')
            ylabel('$\Theta_{t} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsThetaNormToTotalModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsThetaNormToTotalModel.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted, station.theta.meanNormToGsClassModel, idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},nb})$')
            ylabel('$\Theta_{t} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormDistWeightedVsThetaNormToTotalModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormDistWeightedVsThetaNormToTotalModel.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberBedD50, station.theta.meanNormToGsClassModel, idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{50},b})$')
            ylabel('$\Theta_{gs} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberBedD50VsThetaNormToTotalModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberBedD50VsThetaNormToTotalModel.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberBedDistWeighted, station.theta.meanNormToGsClassModel, idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},b})$')
            ylabel('$\Theta_{gs} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberBedDistWeightedVsThetaNormToTotalModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberBedDistWeightedVsThetaNormToTotalModel.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.ustarCalib, station.theta.meanNormToGsClassModel, idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('shear velocity $(u_*)$')
            ylabel('$\Theta_{gs} / c_b$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/UstarVsThetaNormToGsClassModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/UstarVsThetaNormToGsClassModel.eps');
                pause(1);
            end
    end
    
    %% concProf plots -- alpha
    if false 
        % plots with alpha on y axis, attempt to explain variation    
        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormD50, station.alpha.byPredD50RouToTotalModelRou, idx.stationYr, idx.stationOut, f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{50},nb})$')
            ylabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,t})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsAlphaByPredD50RouToTotalModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsAlphaByPredD50RouToTotalModelRou.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted, station.alpha.byPredBulkRouToTotalModelRou, idx.stationYr, idx.stationOut, f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},nb})$')
            ylabel('$\alpha$ from $(Z_{R,p,w} / Z_{R,m,t})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormDistWeightedVsAlphaByPredBulkRouToTotalModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormDistWeightedVsAlphaByPredBulkRouToTotalModelRou.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormD50, station.alpha.byPredD50RouToGsClassModelRou, idx.stationYr, idx.stationOut, f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{50},nb})$')
            ylabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsAlphaByPredD50RouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormD50VsAlphaByPredD50RouToGsClassModelRou.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted, station.alpha.byPredBulkRouToGsClassModelRou, idx.stationYr, idx.stationOut, f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},nb})$')
            ylabel('$\alpha$ from $(Z_{R,p,w} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedNWnormDistWeightedVsAlphaByPredBulkRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedNWnormDistWeightedVsAlphaByPredBulkRouToGsClassModelRou.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberBedD50,  station.alpha.byPredD50RouToTotalModelRou, idx.stationYr, idx.stationOut, f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{50},b})$')
            ylabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,t})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberBedD50VsAlphaByPredD50RouToTotalModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberBedD50VsAlphaByPredD50RouToTotalModelRou.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberBedDistWeighted,  station.alpha.byPredD50RouToTotalModelRou, idx.stationYr, idx.stationOut, f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},b})$')
            ylabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,t})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberBedDistWeightedVsAlphaByPredD50RouToTotalModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberBedDistWeightedVsAlphaByPredD50RouToTotalModelRou.eps');
                pause(1);
            end   

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberBedD50,  station.alpha.byPredD50RouToGsClassModelRou, idx.stationYr, idx.stationOut, f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{50},b})$')
            ylabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberBedD50VsAlphaByPredD50RouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberBedD50VsAlphaByPredD50RouToGsClassModelRou.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberBedDistWeighted,  station.alpha.byPredD50RouToGsClassModelRou, idx.stationYr, idx.stationOut, f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},b})$')
            ylabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberBedDistWeightedVsAlphaByPredD50RouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberBedDistWeightedVsAlphaByPredD50RouToGsClassModelRou.eps');
                pause(1);
            end   

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberBedD50,  station.alpha.byPredBulkRouToGsClassModelRou, idx.stationYr, idx.stationOut, f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{50},b})$')
            ylabel('$\alpha$ from $(Z_{R,p,w} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberBedD50VsAlphaByPredBulkRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberBedD50VsAlphaByPredBulkRouToGsClassModelRou.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberBedDistWeighted,  station.alpha.byPredBulkRouToGsClassModelRou, idx.stationYr, idx.stationOut, f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},b})$')
            ylabel('$\alpha$ from $(Z_{R,p,w} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberBedDistWeightedVsAlphaByPredBulkRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberBedDistWeightedVsAlphaByPredBulkRouToGsClassModelRou.eps');
                pause(1);
            end   

        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(station.cb5S.byTotalModel, station.alpha.byPredD50RouToTotalModelRou, idx.stationYr, idx.stationOut, f);
            plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            ylim([0 2])
            xlabel('$c_{b,t} / S$')
            ylabel('$\alpha (Z_{R,p,50} / Z_{R,m,t})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5STotalModelvsAlphaByPredD50RouToTotalModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5STotalModelvsAlphaByPredD50RouToTotalModelRou.eps');
                pause(1);
            end
            % do a regression to this dataset
            alphaRegress = station.alpha.byPredD50RouToTotalModelRou;
            cbSRegress = station.cb5S.byTotalModel;
            alphaRegressTable = array2table([alphaRegress, cbSRegress], 'VariableNames',  {'alpha', 'cbS'});
            c = 0.9;
            alphaModel = fitnlm(alphaRegressTable, ['alpha ~ 1-(  (a * (cbS) ^ b) / (1 + ((a * (cbS) ^ b) / ',num2str(c), '))  ) '], [0.136, 0.7]);
            plot(1:5:200, 1-(((alphaModel.Coefficients.Estimate(1) .* (1:5:200) .^ alphaModel.Coefficients.Estimate(2)))  ./ ...
                             (1 + (alphaModel.Coefficients.Estimate(1) .* (1:5:200) .^ alphaModel.Coefficients.Estimate(2)) ./ c)), ...
                             'kx-', 'MarkerSize', 4);

        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(station.cb5S.byTotalModel, station.alpha.byPredBulkRouToTotalModelRou, idx.stationYr, idx.stationOut, f);
            plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            ylim([0 2])
            xlabel('$c_{b,t} / S$')
            ylabel('$\alpha (Z_{R,p,w} / Z_{R,m,t})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5STotalModelvsAlphaByBulkD50RouToTotalModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5STotalModelvsAlphaByPredBulkRouToTotalModelRou.eps');
                pause(1);
            end

        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(station.cb5S.byTotalModel, station.alpha.byPredD50RouToGsClassModelRou, idx.stationYr, idx.stationOut, f);
            plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            ylim([0 2])
            xlabel('$c_{b,t} / S$')
            ylabel('$\alpha (Z_{R,p,50} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5STotalModelvsAlphaByPredD50RouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5STotalModelvsAlphaByPredD50RouToGsClassModelRou.eps');
                pause(1);
            end

        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(station.cb5S.byTotalModel, station.alpha.byPredBulkRouToGsClassModelRou, idx.stationYr, idx.stationOut, f);
            plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            ylim([0 2])
            xlabel('$c_{b,t} / S$')
            ylabel('$\alpha (Z_{R,p,w} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5STotalModelvsAlphaByPredBulkRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5STotalModelvsAlphaByPredBulkRouToGsClassModelRou.eps');
                pause(1);
            end

        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(station.cb5S.byGsClassSummedModel, station.alpha.byPredBulkRouToGsClassModelRou, idx.stationYr, idx.stationOut, f);
            plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            ylim([0 2])
            xlabel('$c_{b,\sum i} / S$')
            ylabel('$\alpha$ from $(Z_{R,p,w} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5SgsClassSummedModelvsAlphaByPredBulkRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5SgsClassSummedModelvsAlphaByPredBulkRouToGsClassModelRou.eps');
                pause(1);
            end

        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(station.cb5S.byGsClassSummedModel, station.alpha.byKredMY, idx.stationYr, idx.stationOut, f);
            plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            ylim([0 2])
            xlabel('$c_{b,\sum i} / S$')
            ylabel('$\alpha \equiv K_{red,WP}$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5SgsClassSummedModelvsAlphaByKred.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5SgsClassSummedModelvsAlphaByKred.eps');
                pause(1);
            end

        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(station.cb5S.byGsClassSummedModel(velProfIdx), station.alpha.byKredVelocityModel, idx.stationYr(velProfIdx), idx.stationOut(velProfIdx), f);
            plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            ylim([0 2])
            xlabel('$c_{b,\sum i} / S$')
            ylabel('$\alpha \equiv K_{red,vel}$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5SgsClassSummedModelvsAlphaByKred.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5SgsClassSummedModelvsAlphaByKred.eps');
                pause(1);
            end

        % alpha vs alpha   ////   alpha vs Kred plots
        f = figure(); hold on;
            plotYRdata(station.alpha.byPredBulkRouToGsClassModelRou, station.alpha.byKredMY, idx.stationYr, idx.stationOut, f);
            plot([0 1], [0 1], 'k:', 'LineWidth', 1.5)
            xlabel('$\alpha$ from $(Z_{R,p,w} / Z_{R,m,w})$')
            ylabel('$\alpha \equiv K_{red}$')
            xlim([0 1])
            ylim([0 1])
            axis square
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/alpha.png');
                print('-depsc', '-painters', '-r300', './figsExport/alpha.eps');
                pause(1);
            end
            
        f = figure(); hold on;
            plotYRdata(station.alpha.byPredD50RouToTotalModelRou, station.alpha.byKredMY, idx.stationYr, idx.stationOut, f);
            plot([0 1], [0 1], 'k:', 'LineWidth', 1.5)
            xlabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,t})$')
            ylabel('$\alpha \equiv K_{red}$')
            xlim([0 1])
            ylim([0 1])
            axis square
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/alpha.png');
                print('-depsc', '-painters', '-r300', './figsExport/alpha.eps');
                pause(1);
            end
            
        f = figure(); hold on;
            plotYRdata(station.alpha.byPredD50RouToGsClassModelRou, station.alpha.byKredMY, idx.stationYr, idx.stationOut, f);
            plot([0 1], [0 1], 'k:', 'LineWidth', 1.5)
            xlabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,w})$')
            ylabel('$\alpha \equiv K_{red}$')
            xlim([0 1])
            ylim([0 1])
            axis square
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/alpha.png');
                print('-depsc', '-painters', '-r300', './figsExport/alpha.eps');
                pause(1);
            end
            
        f = figure(); hold on;
            plotYRdata(station.alpha.byPredD50RouToGsClassModelRou, station.alpha.byKredMY, idx.stationYr, idx.stationOut, f);
            plot([0 1], [0 1], 'k-', 'LineWidth', 1.5)
            xlabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,w})$')
            ylabel('$\alpha \equiv K_{red,WP}$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/KredFromWP04ModelVsAlphaByPredD50RouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/KredFromWP04ModelVsAlphaByPredD50RouToGsClassModelRou.eps');
                pause(1);
            end
            
        f = figure(); hold on;
            plotYRdata(station.alpha.byPredD50RouToGsClassModelRou, station.alpha.byKredMY, idx.stationYr, idx.stationOut, f);
            plot([0 1], [0 1], 'k-', 'LineWidth', 1.5)
            xlabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,w})$')
            ylabel('$\alpha \equiv K_{red,WP}$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/KredFromWP04ModelVsAlphaByPredD50RouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/KredFromWP04ModelVsAlphaByPredD50RouToGsClassModelRou.eps');
                pause(1);
            end
            
        f = figure(); hold on;
            plotYRdata(station.alpha.byKredVelocityModel, station.alpha.byKredMY(velProfIdx), idx.stationYr(velProfIdx), idx.stationOut(velProfIdx), f);
            plot([0 2], [0 2], 'k:', 'LineWidth', 1.5)
            xlabel('$\alpha \equiv \alpha_{VEL}$')
            ylabel('$\alpha \equiv K_{red}$')
            axis square
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/alpha.png');
                print('-depsc', '-painters', '-r300', './figsExport/alpha.eps');
                pause(1);
            end
            
        f = figure(); hold on;
            plotYRdata(station.alpha.byKredVelocityModel, station.alpha.byPredBulkRouToGsClassModelRou(velProfIdx), idx.stationYr(velProfIdx), idx.stationOut(velProfIdx), f);
            plot([0 2], [0 2], 'k:', 'LineWidth', 1.5)
            xlabel('$\alpha \equiv \alpha_{VEL}$')
            ylabel('$\alpha \equiv \alpha_f$')
            axis square
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/alpha.png');
                print('-depsc', '-painters', '-r300', './figsExport/alpha.eps');
                pause(1);
            end
            
        f = figure(); hold on;
            plotYRdata(station.alpha.byWP04func, station.alpha.byKredMY, idx.stationYr, idx.stationOut, f);
            plot([0 2], [0 2], 'k:', 'LineWidth', 1.5)
            xlabel('$\alpha \equiv \alpha_{WP04}$')
            ylabel('$\alpha \equiv K_{red}$')
            axis square
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/alpha.png');
                print('-depsc', '-painters', '-r300', './figsExport/alpha.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.alpha.byPredBulkRouToGsClassModelRou, station.alpha.byPredBulkRouToTotalModelRou, idx.stationYr, idx.stationOut, f);
            plot([0 2], [0 2], 'k:', 'LineWidth', 1.5)
            xlabel('$\alpha$ from $(Z_{R,p,w} / Z_{R,m,w})$')
            ylabel('$\alpha$ from $(Z_{R,p,w} / Z_{R,m,t})$')
            axis square
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/.png');
                print('-depsc', '-painters', '-r300', './figsExport/.eps');
                pause(1);
            end

        f = figure(); hold on; %#ok<*NASGU>
            plot(station.alpha.byPredBulkRouToGsClassModelRou(and(velProfIdx, station.year==2016)), station.alpha.byKredMY(station.velCalib.ustarYr==2016), ...
                'Marker', mar2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col2, 'MarkerSize', 6, 'LineStyle', 'none')
            plot(station.alpha.byPredBulkRouToGsClassModelRou(and(velProfIdx, station.year==2018)), station.alpha.byKredMY(station.velCalib.ustarYr==2018), ... 
                'Marker', mar3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col3, 'MarkerSize', 6, 'LineStyle', 'none')
            plot([0 2], [0 2], 'k-', 'LineWidth', 1);
            axis square
            xlabel('$\alpha$ from $(Z_{R,p,w} / Z_{R,m,w})$')
            ylabel('$\alpha \equiv K_{red,WP}$')
            box on
            legend('2016', '2018', 'Location', 'NorthWest')
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/alphaByKredVsalphaByConcProf.png');
                print('-depsc', '-painters', '-r300', './figsExport/alphaByKredVsalphaByConcProf.eps');
                pause(1);
            end

        % plot the rouse numbers vs one another to see whats up with the differences
        f = figure(); hold on;
            subplot(2, 2, 1)
                plotYRdata(station.Rou.gsClassModelBulkRou, station.Rou.gsClassPredBulkRou, idx.stationYr, idx.stationOut, f);
                plot([0 1], [0 1], 'k:', 'LineWidth', 1.5)
                xlabel('$(Z_{R,m,w})$')
                ylabel('$(Z_{R,p,w})$')
                axis square
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            subplot(2, 2, 2)
                plotYRdata(station.Rou.gsClassModelBulkRou, station.Rou.gsClassPredD50Rou, idx.stationYr, idx.stationOut, f);
                plot([0 1], [0 1], 'k:', 'LineWidth', 1.5)
                xlabel('$(Z_{R,m,w})$')
                ylabel('$(Z_{R,p,50})$')
                axis square
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            subplot(2, 2, 3)
                plotYRdata(station.Rou.gsClassModelBulkRou, station.Rou.totalModelRou, idx.stationYr, idx.stationOut, f);
                plot([0 1], [0 1], 'k:', 'LineWidth', 1.5)
                xlabel('$(Z_{R,m,w})$')
                ylabel('$(Z_{R,m,t})$')
                axis square
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            subplot(2, 2, 4)
                plotYRdata(station.Rou.totalModelRou, station.Rou.gsClassPredD50Rou, idx.stationYr, idx.stationOut, f);
                plot([0 1], [0 1], 'k:', 'LineWidth', 1.5)
                xlabel('$(Z_{R,m,t})$')
                ylabel('$(Z_{R,p,50})$')
                axis square
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 850 600], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/RouseNumbersVsRouseNumbers.png');
                print('-depsc', '-painters', '-r300', './figsExport/RouseNumbersVsRouseNumbers.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted .* station.suspension.nearBedNWnormDistWeightedRep .^ 0.6, station.alpha.byPredBulkRouToGsClassModelRou, ...
                idx.stationYr, idx.stationOut, f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity $\times$ Rep $([u_* / w_{s,D_{w},nb}] \cdot $Rep$^{0.6})$')
            ylabel('$\alpha$ from $(Z_{R,p,w} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedDistWeightedRepVsAlphaByPredBulkRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedDistWeightedRepVsAlphaByPredBulkRouToGsClassModelRou.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted, station.Risr.byBedD50, idx.stationYr, idx.stationOut, f);
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},nb}$')
            ylabel('$Ri_{sr,50,b}$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedDistWeightedVsRisrByBedD50.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedDistWeightedVsRisrByBedD50.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted, station.Risr.byBedBulk, idx.stationYr, idx.stationOut, f);
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},nb}$')
            ylabel('$Ri_{sr,w,b}$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedDistWeightedVsRisrByBedBulk.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedDistWeightedVsRisrByBedBulk.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberNearBedNWnormDistWeighted, station.Risr.byNearBedBulk, idx.stationYr, idx.stationOut, f);
            xlabel('dimensionless shear velocity $(u_* / w_{s,D_{w},nb}$')
            ylabel('$Ri_{sr,w,nb}$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedDistWeightedVsRisrByNearBedBulk.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedDistWeightedVsRisrByNearBedBulk.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.Risr.byBedBulk, station.alpha.byPredD50RouToGsClassModelRou, idx.stationYr, idx.stationOut, f);
            xlabel('$Ri_{sr,w,b}$')
            ylabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,w})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/RisrByBedBulkVsAlphaByPredD50RouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/RisrByBedBulkVsAlphaByPredD50RouToGsClassModelRou.eps');
                pause(1);
            end

        f = figure(); hold on;
            plotYRdata(station.alpha.byPredD50RouToGsClassModelRou, station.Risr.byNearBedBulk, idx.stationYr, idx.stationOut, f);
            xlabel('$\alpha$ from $(Z_{R,p,50} / Z_{R,m,w})$')
            ylabel('$Ri_{sr,w,nb}$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/RisrByNearBedBulkVsAlphaByPredD50RouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/RisrByNearBedBulkVsAlphaByPredD50RouToGsClassModelRou.eps');
                pause(1);
            end
    end
       
        
    %% plots with theta/alpha on y axis, for the grain size specific profile fits
    if false
        nClasses = size(df(1).nearBedData.gsDistNearBedNWnorm, 1);
        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberGsClass, station.theta.meanNormGsClassToGsClassModel, ...
                varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f);
            plot(xlim, [0 0], 'k:', 'LineWidth', 1.5)
            xlabel('dimensionless shear velocity ($u_* / w_{s,i}$)')
            ylabel('$\Theta_i / c_{b,i}$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedGsClassVsThetaMeanNormGsClassToGsClassModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedGsClassVsThetaMeanNormGsClassToGsClassModel.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberGsClass, station.alpha.byPredGsClassRouToGsClassModelRou, ...
                varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f);
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            text(0.05, 0.9, '\^outliers', 'Units', 'normalized')
            ylim([0 5])
            xlabel('dimensionless shear velocity $(u_* / w_{s,i})$')
            ylabel('$\alpha$ from $(Z_{R,p,i} / Z_{R,m,i})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedGsClassVsAlphaByPredGsClassRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedGsClassVsAlphaByPredGsClassRouToGsClassModelRou.eps');
                pause(1);
            end    

        f = figure(); hold on;
            plotYRdata(station.suspension.suspensionNumberGsClass, station.cb5S.byGsClassModel, ...
                varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f);
            xlabel('dimensionless shear velocity $(u_* / w_{s,i})$')
            ylabel('$c_{b,i} / S$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/suspensionNumberNearBedGsClassVscb5SByGsClassModel.png');
                print('-depsc', '-painters', '-r300', './figsExport/suspensionNumberNearBedGsClassVscb5SByGsClassModel.eps');
                pause(1);
            end    

        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(station.cb5S.byGsClassModel, station.alpha.byPredGsClassRouToGsClassModelRou, ...
                varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f);
            plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            text(0.05, 0.9, '\^outliers', 'Units', 'normalized')
            ylim([0 2])
            xlabel('$c_{b,i} / S$')
            ylabel('$\alpha$ from $(Z_{R,p,i} / Z_{R,m,i})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5SGsClasslModelvsAlphaByPredGsClassRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5SGsClasslModelvsAlphaByPredGsClassRouToGsClassModelRou.eps');
                pause(1);
            end
            
        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(station.cb5S.byGsClassModelToGsClassFrac, station.alpha.byPredGsClassRouToGsClassModelRou, ...
                varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f);
            plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            text(0.05, 0.9, '\^outliers', 'Units', 'normalized')
            ylim([0 2])
            xlabel('$c_{b} / S$')
            ylabel('$\alpha$ from $(Z_{R,p,i} / Z_{R,m,i})$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5SGsClasslModelvsAlphaByPredGsClassRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5SGsClasslModelvsAlphaByPredGsClassRouToGsClassModelRou.eps');
                pause(1);
            end
            
            
        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(station.suspension.suspensionNumberGsClass, station.alpha.byPredGsClassRouToGsClassModelRou, ...
                varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, 'ColorByVariable', concatenate_cell_vectors(station.suspension.nearBedNWnormDist));
            % plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            % text(0.05, 0.9, '\^{}outliers', 'Units', 'normalized')
            ylim([0 1000])
            xlabel('$u*/ws$')
            ylabel('$\alpha$ from $(Z_{R,p,i} / Z_{R,m,i})$')
            set(gca, 'xscale', 'log', 'yscale', 'log')
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            colorbar;
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5SGsClasslModelvsAlphaByPredGsClassRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5SGsClasslModelvsAlphaByPredGsClassRouToGsClassModelRou.eps');
                pause(1);
            end
            
            
        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata(concatenate_cell_vectors(station.suspension.nearBedNWnormDistVals)./100.*...
                (1./station.suspension.suspensionNumberGsClass).*(station.cb5S.byGsClassModelToGsClassFrac), ...
                station.alpha.byPredGsClassRouToGsClassModelRou, ...
                varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, 'ColorByVariable', concatenate_cell_vectors(station.suspension.nearBedNWnormDist));
            % plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            % text(0.05, 0.9, '\^{}outliers', 'Units', 'normalized')
            ylim([0.01 1000])
            xlim([1, 1e6])
            xlabel('$(u_*/w_s)(cb/S)$')
            ylabel('$\alpha$ from $(Z_{R,p,i} / Z_{R,m,i})$')
            set(gca, 'xscale', 'log', 'yscale', 'log')
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            colorbar;
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5SGsClasslModelvsAlphaByPredGsClassRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5SGsClasslModelvsAlphaByPredGsClassRouToGsClassModelRou.eps');
                pause(1);
            end
            
        f = figure(); hold on;
            fill(mods.WPdatacloud(:,1), mods.WPdatacloud(:,2), [0.8, 0.8, 0.8], 'EdgeColor', 'none');
            plotYRdata((station.cb5S.byGsClassModel), station.alpha.byPredGsClassRouToGsClassModelRou, ...
                varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, 'ColorByVariable', concatenate_cell_vectors(station.suspension.nearBedNWnormDist));
            plot(mods.WPalphaEvalXs, mods.WPalphaEval, 'k-x', 'MarkerSize', 2)
            % text(0.05, 0.9, '\^{}outliers', 'Units', 'normalized')
            % ylim([0.01 1000])
            % xlim([1, 1e6])
            xlabel('$c_{b,i} / S$')
            ylabel('$\alpha$ from $(Z_{R,p,i} / Z_{R,m,i})$')
            set(gca, 'xscale', 'log', 'yscale', 'log')
            plot(xlim, [1 1], 'k:', 'LineWidth', 1.5)
            colorbar;
            set(gca, 'colorscale', 'log')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 250], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
            if printOut
                print('-dpng', '-r300', './figsExport/cb5SGsClasslModelvsAlphaByPredGsClassRouToGsClassModelRou.png');
                print('-depsc', '-painters', '-r300', './figsExport/cb5SGsClasslModelvsAlphaByPredGsClassRouToGsClassModelRou.eps');
                pause(1);
            end
            
        
    end
    
    
    %% every concProf
    if false
        [sortByUstar, sortByUstarIdx] = sort( cell2mat(arrayfun(@(s) s.Velocity.ustarCalib, df, 'Unif', 0)) ); %#ok<*ASGLU>
        plotMatLen = ceil(sqrt(size(df,1)));
        m_plot = 0.025; % margins
        s_plot = 0.0125; % spacing
        w_plot = (1-2*m_plot-6*s_plot) / 7;
        h_plot = (1-2*m_plot-7*s_plot) / 8;
        x_pos = repmat(m_plot+(0:6)*w_plot+(0:6)*s_plot, 8, 1)';
        y_pos = repmat(1-(h_plot+m_plot+(0:7)*h_plot+(0:7)*s_plot), 7, 1);
        figure()
        for p = 1:size(df,1)
            dfIdx = sortByUstarIdx(p);
            if datenum(df(dfIdx).CollectionDate) < datenum('01/01/2016') 
                yrCol = col1;
                yrMark = mar1;
            elseif datenum(df(dfIdx).CollectionDate) < datenum('01/01/2017')
                yrCol = col2;
                yrMark = mar2;
            else
                yrCol = col3;
                yrMark = mar3;
            end           
            subplot('Position', [x_pos(p) y_pos(p) w_plot, h_plot]); hold on;
                plot(df(dfIdx).waterSamplesTable.concNW, df(dfIdx).waterSamplesTable.sampleZ, 'k', 'MarkerFaceColor', yrCol, 'MarkerSize', 4, 'Marker', yrMark, 'LineStyle', 'none')
%                 plot(df(dfIdx).concProf.totalModel.Cs, df(dfIdx).concProf.totalModel.Zs, 'k--', 'MarkerSize', 3)
                plot(df(dfIdx).concProf.gsClassModel.CsSum, df(dfIdx).concProf.gsClassModel.ZsSum, 'k--', 'MarkerSize', 3)
                plot(df(dfIdx).concProf.gsClassPred.CsSum, df(dfIdx).concProf.gsClassPred.ZsSum, 'k-', 'MarkerSize', 3)
                plot(df(dfIdx).concProf.WPalphaPred.CsSum, df(dfIdx).concProf.WPalphaPred.ZsSum, 'b-', 'MarkerSize', 3)
                MYsum_mdl = plot(df(dfIdx).concProf.MYfullPred.CsSum, df(dfIdx).concProf.MYfullPred.Zs, ...
                    'r--', 'LineWidth', 1.5); % wright parker prediction
                box on; set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
                xl = [0 max(df(dfIdx).waterSamplesTable.concNW)*1.3];
                xlim(xl);
                yl = [0 df(dfIdx).FlowDepthAtCollection*1.2];
                ylim(yl);
                set(gca, 'xTick', (0:floor(xl(end))))
                set(gca, 'yTick', (0:0.5:floor(yl(end))))
                xText = df(dfIdx).concProf.gsClassPred.CsSum( floor(length(df(dfIdx).concProf.totalModel.Cs)/2) )*1.05;
                text(xText, yl(end)*0.475, ['$c_b$ = ', num2str(round(df(dfIdx).concProf.totalModel.Cs(1),2)), ' g/L'])
                text(xText, yl(end)*0.875, df(dfIdx).StationID, 'Interpreter', 'none', 'FontSize', 10)
                text(xText, yl(end)*0.675, ['$u_*$ = ', num2str(round(df(dfIdx).Velocity.ustarCalib*100,2)), ' cm/s'])
                text(xl(end)*0.1, yl(end)*0.1, ['H = ', num2str(round(df(dfIdx).FlowDepthAtCollection,1)), ' m'], 'rotation', 90)
                set(gca, 'xticklabels', [])
                set(gca, 'yticklabels', [])
                set(gca,'TickLength',[0.015, 0.01])
        end
        p = p + 1;
        subplot('Position', [x_pos(p) y_pos(p) w_plot, h_plot])
            xlabel('conc (g/L)')
            ylabel('dist. above bed (m)')
            legend('samples', 'best-fit', 'pred.')
            box on; set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gca, 'xticklabels', [])
            set(gca, 'yticklabels', [])
            set(gca,'TickLength',[0.015, 0.01])
        set(gcf, 'Pos', [50 100 1600 1500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
        if printOut
            print('-dpng', '-r300', './figsExport/everyConcProf.png');
            print('-depsc', '-painters', '-r300', './figsExport/everyConcProf.eps');
            pause(1);
        end
    end
    
    % alpha regressions for grain size bins
    if false
        cmap = parula(nClasses); %flipud(jet(6));
        m_plot = 0.1; % margins
        s_plot = 0.125; % spacing
        nx = 2;
        ny = 2;
        w_plot = (1-2*m_plot-(nx-1)*s_plot) / nx;
        h_plot = (1-2*m_plot-(ny-1)*s_plot) / ny;
        x_pos = repmat(m_plot+(0:(nx-1))*w_plot+(0:(nx-1))*s_plot, ny, 1)';
        y_pos = repmat(1-(h_plot+m_plot+(0:(ny-1))*h_plot+(0:(ny-1))*s_plot), nx, 1);

        figure();
        subplot('position', [x_pos(1) y_pos(1) w_plot, h_plot]); hold on
            plot([0 7], [0 0], '-', 'Color', [0.6 0.6 0.6])
            boxplot(squeeze(reg.msd)')
            cleanup_boxplot('color', cmap)
            set(gca, 'xticklabels', arrayfun(@(x) num2str(round(x,0)), gsDistClass2num(dfd(1).nearBedData.gsDistNearBedNWnorm), 'unif', 0))
            xlabel('grain-size ($\mu$m)')
            ylabel('msd')
        subplot('position', [x_pos(3) y_pos(3) w_plot, h_plot]); hold on
            plot([0 7], [0 0], '-', 'Color', [0.6 0.6 0.6])
            boxplot(squeeze(reg.nmsd)')
            cleanup_boxplot('color', cmap)
            set(gca, 'xticklabels', arrayfun(@(x) num2str(round(x,0)), gsDistClass2num(dfd(1).nearBedData.gsDistNearBedNWnorm), 'unif', 0))
            xlabel('grain-size ($\mu$m)')
            ylabel('nmsd')
        subplot('position', [x_pos(2) y_pos(2) w_plot, h_plot]); hold on
            plot([0 max(reg.Zs)], [0 0], '-', 'Color', [0.6 0.6 0.6])
            sddstd_h = fill(repmat([reg.Zs(1:49); flipud(reg.Zs(1:49))], 1, 3), [reg.sdd(1:49,1:3)+reg.sddstd(1:49,1:3); flipud(reg.sdd(1:49,1:3)-reg.sddstd(1:49,1:3))], [0 0 0], 'LineWidth', 1);
            sdd_h = plot(reg.Zs(1:49), reg.sdd(1:49,1:3), 'LineWidth', 2);
            set(sdd_h, {'color'}, num2cell(cmap(1:3,:), 2));
            set(sddstd_h, {'facecolor'}, num2cell(cmap(1:3,:), 2), {'edgecolor'}, num2cell(cmap(1:3,:), 2), 'facealpha', 0.2);
            xlim([0 max(reg.Zs)])
            xlabel('z/H (-)')
            cleanup_boxplot('color', jet(0))
        subplot('position', [x_pos(4) y_pos(4) w_plot, h_plot]); hold on
            plot([0 max(reg.Zs)], [0 0], '-', 'Color', [0.6 0.6 0.6])
            sddstd_h = fill(repmat([reg.Zs(1:49); flipud(reg.Zs(1:49))], 1, 3), [reg.sdd(1:49,4:6)+reg.sddstd(1:49,4:6); flipud(reg.sdd(1:49,4:6)-reg.sddstd(1:49,4:6))], [0 0 0], 'LineWidth', 1);
            sdd_h = plot(reg.Zs(1:49), reg.sdd(1:49,4:6), 'LineWidth', 2);
            set(sdd_h, {'color'}, num2cell(cmap(4:6,:), 2));
            set(sddstd_h, {'facecolor'}, num2cell(cmap(4:6,:), 2), {'edgecolor'}, num2cell(cmap(4:6,:), 2), 'facealpha', 0.2);
            xlim([0 max(reg.Zs)])
            xlabel('z/H (-)')
            cleanup_boxplot('color', jet(0))


         figure();
         subplot(2, 2, [1 3]); hold on
            plot([0 0], [0 max(reg.Zs)], '-', 'Color', [0.6 0.6 0.6])
            % sddstd_h = fill([sdd(1:49,:)+sddstd(1:49,:); flipud(sdd(1:49,:)-sddstd(1:49,:))].*2650, repmat([Zs(1:49); flipud(Zs(1:49))], 1, 6), [0 0 0], 'LineWidth', 0.5);
            sdd_h = plot(reg.sdd(1:49,:).*2650, reg.Zs(1:49), 'LineWidth', 2);
            set(sdd_h, {'color'}, num2cell(cmap, 2));
            % set(sddstd_h, {'facecolor'}, num2cell(cmap, 2), {'edgecolor'}, num2cell(cmap, 2), 'facealpha', 0.1);
            ylim([0 max(reg.Zs)])
            xlim([-0.5 0.5])
            ylabel('z/H (-)')
            xlabel('msd (g/L)')
            cleanup_boxplot('color', jet(0))
         subplot('position', [x_pos(2) y_pos(2) w_plot, h_plot]); hold on
            plot([0 7], [0 0], '-', 'Color', [0.6 0.6 0.6])
            boxplot(squeeze(reg.msd)'.*2650)
            cleanup_boxplot('color', cmap)
            set(gca, 'xticklabels', arrayfun(@(x) num2str(round(x,0)), gsDistClass2num(dfd(1).nearBedData.gsDistNearBedNWnorm), 'unif', 0))
            xlabel('grain-size ($\mu$m)')
            ylabel('msd (g/L)')
        subplot('position', [x_pos(4) y_pos(4) w_plot, h_plot]); hold on
            plot([0 7], [0 0], '-', 'Color', [0.6 0.6 0.6])
            boxplot(squeeze(reg.nmsd)')
            cleanup_boxplot('color', cmap)
            set(gca, 'xticklabels', arrayfun(@(x) num2str(round(x,0)), gsDistClass2num(dfd(1).nearBedData.gsDistNearBedNWnorm), 'unif', 0))
            xlabel('grain-size ($\mu$m)')
            ylabel('nmsd (-)')
        
    end
    
    %% PROFILE COMPARIONS
    if false
        %% theta by grainsize plots: MYfull - datafits
        cmap = parula(nClasses); 
        m_plot = 0.1; % margins
        s_plot = 0.125; % spacing
        nx = 2;
        ny = 2;
        w_plot = (1-2*m_plot-(nx-1)*s_plot) / nx;
        h_plot = (1-2*m_plot-(ny-1)*s_plot) / ny;
        x_pos = repmat(m_plot+(0:(nx-1))*w_plot+(0:(nx-1))*s_plot, ny, 1)';
        y_pos = repmat(1-(h_plot+m_plot+(0:(ny-1))*h_plot+(0:(ny-1))*s_plot), nx, 1);
        figure();
        subplot(2, 2, [1 3]); hold on
            plot([0 0], [0 max(sdif.thet.Zs)], '-', 'Color', [0.6 0.6 0.6])
            % sddstd_h = fill([sdd(1:49,:)+sddstd(1:49,:); flipud(sdd(1:49,:)-sddstd(1:49,:))].*2650, repmat([Zs(1:49); flipud(Zs(1:49))], 1, 6), [0 0 0], 'LineWidth', 0.5);
            sdd_h = plot(sdif.thet.sdd(1:49,:).*2650, sdif.thet.Zs(1:49), 'LineWidth', 2);
            set(sdd_h, {'color'}, num2cell(cmap, 2));
            % set(sddstd_h, {'facecolor'}, num2cell(cmap, 2), {'edgecolor'}, num2cell(cmap, 2), 'facealpha', 0.1);
            text(0.1, 0.95, 'a', 'units', 'normalized', 'FontSize', 14)
            ylim([0 max(sdif.thet.Zs)])
            xlim([-0.5 0.5])
            ylabel('$z/H$ (-)')
            xlabel('$\xi_{WP04,i}$ (g/L)')
            cleanup_boxplot('color', jet(0))
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        subplot('position', [x_pos(2) y_pos(2) w_plot, h_plot]); hold on
            plot([0 7], [0 0], '-', 'Color', [0.6 0.6 0.6])
            text(0.9, 0.9, 'b', 'units', 'normalized', 'FontSize', 14)
            boxplot(squeeze(sdif.thet.msd)'.*2650)
            cleanup_boxplot('color', cmap)
            set(gca, 'xticklabels', arrayfun(@(x) num2str(round(x,0)), gsDistClass2num(df(1).nearBedData.gsDistNearBedNWnorm), 'unif', 0))
            xlabel('grain-size ($\mu$m)')
            ylabel('$\int_b^H \xi_{WP04,i}/H$ (g/L)')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        subplot('position', [x_pos(4) y_pos(4) w_plot, h_plot]); hold on
            plot([0 7], [0 0], '-', 'Color', [0.6 0.6 0.6])
            text(0.9, 0.9, 'c', 'units', 'normalized', 'FontSize', 14)
            boxplot(squeeze(sdif.thet.nmsd)')
            cleanup_boxplot('color', cmap)
            set(gca, 'xticklabels', arrayfun(@(x) num2str(round(x,0)), gsDistClass2num(df(1).nearBedData.gsDistNearBedNWnorm), 'unif', 0))
            xlabel('grain-size ($\mu$m)')
            ylabel('$[\int_b^H \xi_{WP04,i}/H]/\bar{c}_b$ (-)')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
            set(gcf, 'Pos', [50 100 600 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
        if printOut
            print('-dpng', '-r300', [print_root, 'WPalpha_grainsize_msd.png']);
            print('-depsc', '-painters', '-r300', [print_root, 'WPalpha_grainsize_msd.eps']);
            pause(1);
        end
    end
    
    
    %% entrainment plots
    % entrainment with WP04 relation by grain size class
    if false
        f = figure();
        hold on; set(gca, 'XScale', 'log', 'YScale', 'log')
            lineGP91 = loglog(entr.mods.GP91Zu, entr.mods.GP91Es, 'Color', [0.7 0.7 0.7]);
            lineWP04 = loglog(entr.mods.WP04Zu, entr.mods.WP04Es, 'Color', [0 0 0]);
            plotYRdata(entr.dist.lXiLong, entr.dist.EsiMeasLong, repelem(idx.sampleYr, 6, 1), repelem(idx.sampleOut, 6, 1), f, ...
                'ColorByVariable', concatenate_cell_vectors(varRepLong(arrayfun(@(x) gsDistClass2num(x.bedData.gsDistBed), df, 'unif', 0), nearBed.dims)))
            text(0.8, 0.1, 'color = D', 'units', 'normalized')
            ylim([1e-7 1e0])
            xlabel(['$\lambda X_i = ( \frac{u_{*,sk}}{w_{si}} Re_{pi}^{0.6} ) S^{' num2str(slope_exp) '} (\frac{D_i}{D_{50}} )^{0.2}$'], 'Interpreter', 'latex')
            ylabel('$E_{si} \equiv \bar{c}_b / F_i$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
            set(gcf, 'Pos', [100 100 400 600]);
            set(gcf, 'PaperPositionMode', 'auto')
            drawnow; pause(0.1);
            if printOut
                print('-dpng', '-r300', './figsExport/entrXiVsEs.png');
                print('-depsc', '-painters', '-r300', './figsExport/entrXiVsEs.eps');
                pause(1);
            end

        f = figure();
        hold on; set(gca, 'XScale', 'log', 'YScale', 'log')
            lineGP91 = loglog(entr.mods.GP91Zu, entr.mods.GP91Es, 'Color', [0.7 0.7 0.7]);
            lineWP04 = loglog(entr.mods.WP04Zu, entr.mods.WP04Es, 'Color', [0 0 0]);
            plotYRdata(entr.dist.lXiLong, entr.dist.EsiMeasLong, repelem(idx.sampleYr, 6, 1), repelem(idx.sampleOut, 6, 1), f, ...
                'ColorByVariable', concatenate_cell_vectors(cellfun(@(ustar, ws) ustar./ws, entr.dist.ustarsk, entr.dist.ws, 'Unif', 0)))
            text(0.8, 0.1, 'color = u*/ws', 'units', 'normalized')
            ylim([1e-7 1e0])
            xlabel(['$\lambda X_i = ( \frac{u_{*,sk}}{w_{si}} Re_{pi}^{0.6} ) S^{' num2str(slope_exp) '} (\frac{D_i}{D_{50}} )^{0.2}$'], 'Interpreter', 'latex')
            ylabel('$E_{si} \equiv \bar{c}_b / F_i$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
            set(gcf, 'Pos', [100 100 400 600]);
            set(gcf, 'PaperPositionMode', 'auto')
            drawnow; pause(0.1);
            if printOut
                print('-dpng', '-r300', './figsExport/entrXiVsEs.png');
                print('-depsc', '-painters', '-r300', './figsExport/entrXiVsEs.eps');
                pause(1);
            end

        f = figure();
        hold on; set(gca, 'XScale', 'log', 'YScale', 'log')
            plotYRdata(entr.dist.ustarskWSRep06long, entr.dist.EsiMeasLong, repelem(idx.sampleYr, 6, 1), repelem(idx.sampleOut, 6, 1), f, ...
                'ColorByVariable', concatenate_cell_vectors(varRepLong(arrayfun(@(x) gsDistClass2num(x.bedData.gsDistBed), df, 'unif', 0), nearBed.dims)))
            text(0.8, 0.1, 'color = D', 'units', 'normalized')
            ylim([1e-7 1e0])
            xlabel(['$ ( \frac{u_{*,sk}}{w_{si}} Re_{pi}^{0.6} ) $'], 'Interpreter', 'latex')
            ylabel('$E_{si} \equiv \bar{c}_b / F_i$')
            box on
            set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
            set(gcf, 'Pos', [100 100 400 600]);
            set(gcf, 'PaperPositionMode', 'auto')
            drawnow; pause(0.1);
            if printOut
                print('-dpng', '-r300', './figsExport/entrustarwsRep06VsEs.png');
                print('-depsc', '-painters', '-r300', './figsExport/entrustarwsRep06VsEs.eps');
                pause(1);
            end

        % entrainment with WP04 relation integrated to each station
    %     f = figure();
    %     hold on; set(gca, 'XScale', 'log', 'YScale', 'log')
    %         lineGP91 = loglog(entr.mods.GP91Zu, entr.mods.GP91Es, 'Color', [0.7 0.7 0.7]);
    %         lineWP04 = loglog(entr.mods.WP04Zu, entr.mods.WP04Es, 'Color', [0 0 0]);
    %         plotYRdata(entr.dist.lXi, entr.dist.EsiMeas, idx.sampleYr, idx.sampleOut, f)
    %         ylim([1e-4 1e0])
    %         xlabel(['$\lambda X_i = ( \frac{u_{*,sk}}{w_{si}} Re_{pi}^{0.6} ) S^{' num2str(slope_exp) '} (\frac{D_i}{D_{50}} )^{0.2}$'], 'Interpreter', 'latex')
    %         ylabel('E_s')
    %         box on
    %         set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
    %         set(gcf, 'Pos', [100 100 400 600]);
    %         set(gcf, 'PaperPositionMode', 'auto')
    %         drawnow; pause(0.1);
    %         if printOut
    %             print('-dpng', '-r300', './figsExport/entrXiVsEs.png');
    %             print('-depsc', '-painters', '-r300', './figsExport/entrXiVsEs.eps');
    %             pause(1);
    %         end

        % WP04 prediction vs measured
    %     f = figure();
    %     hold on; set(gca, 'XScale', 'log', 'YScale', 'log')
    %         plotYRdata(nearBed.concNWVol, entr.dist.EsipredEff, idx.sampleYr, idx.sampleOut, f)
    %         plot([1e-4 1e0], [1e-4 1e0], 'k-', 'LineWidth', 1);
    %         axis square
    %         xlabel('measured conc')
    %         ylabel('pred conc')
    %         box on; grid on;
    %         set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k')
    %         set(gcf, 'Pos', [100 100 400 400]);
    %         set(gcf, 'PaperPositionMode', 'auto')
    %         drawnow; pause(0.1);
    %         if printOut
    %             print('-dpng', '-r300', './figsExport/entrEsMeasVsEsPredWP04.png');
    %             print('-depsc', '-painters', '-r300', './figsExport/entrEsMeasVsEsPredWP04.eps');
    %             pause(1);
    %         end
    end
    
    %% r0 and transport explor
    if false
        f = figure();
        plotYRdata(station.alpha.byPredGsClassRouToGsClassModelRou, station.beta.fieldr0Long, ...
            varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, ...
            'ColorByCategorical', station.gsClassVectorCats)
        set(gca, 'xscale', 'log', 'yscale', 'log')
        xlabel('$alpha_f$ or $K_{red}$')
        ylabel('r0')
        
        
        f = figure();
        plotYRdata(station.alpha.byPredBulkRouToGsClassModelRou(idx.velProf), station.transport.fieldNoStrat_ratio, ...
            idx.stationYr(idx.velProf), idx.stationOut(idx.velProf), f)
       
        
        f = figure();
        plotYRdata(station.transport.pred_Rou, station.beta.fieldr0Long, ...
            varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, ...
            'ColorByCategorical', station.gsClassVectorCats)
            set(gca, 'xscale', 'log')
        % theroretical model
        Rous = logspace(-1, 1, 100);
        Zs = linspace(0.05, 1, 50)';
        Profs = 1 .* ( ((1-Zs)./Zs) ./ 19 ) .^ Rous;
        Cs = mean(Profs,1);
        r0s = 1./Cs;
        
            
                    
        f = figure();
        plotYRdata((station.suspension.nearBedtstarGsClass), station.beta.fieldr0Long.*repelem(station.Cf, nClasses), ...
            varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, ...
            'ColorByCategorical', station.gsClassVectorCats)
            set(gca, 'xscale', 'log', 'yscale', 'log')
    
        f = figure();
        plotYRdata((station.suspension.suspensionNumberUstarskGsClass), station.beta.fieldr0Long.*repelem(station.Cf, nClasses), ...
            varRepLong(idx.stationYr, nClasses), varRepLong(idx.stationOut, nClasses), f, ...
            'ColorByCategorical', station.gsClassVectorCats)
            set(gca, 'xscale', 'log', 'yscale', 'log')
            
%             qs = qs_n .* sqrt(con.Rr * con.g * (D50*D50*D50))

    f = figure();
    hold on; set(gca, 'XScale', 'log', 'YScale', 'log')
        plot([1e-5, 100], [1e-5, 100], 'k-')
        plot([1e-5, 100], [1e-6, 10], 'k--')
        plot([1e-6, 100], [1e-5, 1000], 'k--')
        plotYRdata(entr.dist.EsiMeasLong, concatenate_cell_vectors(varRepLong(num2cell(reshape(station.transport.Esolve, [], nClasses),2), nearBed.dims)), repelem(idx.sampleYr, 6, 1), repelem(idx.sampleOut, 6, 1), f, ...
            'ColorByVariable', concatenate_cell_vectors(varRepLong(arrayfun(@(x) repmat(x, 1, nClasses), ...
            station.alpha.byPredBulkRouToGsClassModelRou, 'Unif', 0), nearBed.dims)), 'IgnoreOutliers', true);

    end

    
    %% richardson(-like) numbers plots
    

    %% 2015 vs 2016 vs 18 plots

    % initialize long versions of table for pulling apart as boxplots
    longWaterTable = []; % long table of all water sample concentrations
    longWaterCollDate = []; % matching length table of all water sample dates

    % set up a station table of extracted variables
    stationTableVars = {'CollectionDate', 'gsWLcutoff', 'gsSummBedD05', 'gsSummBedD50', 'gsSummBedD90'};
    stationTable = cell2table(cell(0, length(stationTableVars)), 'VariableNames', stationTableVars);

    % loop through datatable to build up long versions of tables (initialized above)
    for i = 1:size(df, 1)
        longWaterTable = [longWaterTable; df(i).waterSamplesTable];
        longWaterCollDate = [longWaterCollDate; repmat(datenum(df(i).CollectionDate), size(df(i).waterSamplesTable,1), 1)];
        stationTable = [stationTable; {df(i).CollectionDate, df(i).Samples(1).gsWLdata.gsWLcutoff,  df(i).bedData.gsSummBed.d5, df(i).bedData.gsSummBed.d50, df(i).bedData.gsSummBed.d90}];
    end

    % separate out the long tables into year by year data
    waterTable15 = longWaterTable(longWaterCollDate < datenum('01/01/2016'), :);
    waterTable16 = longWaterTable(and(longWaterCollDate > datenum('01/01/2016'), longWaterCollDate < datenum('01/01/2017')), :);
    waterTable18 = longWaterTable(longWaterCollDate > datenum('01/01/2018'), :);

    bedTable15 = stationTable(datenum(stationTable.CollectionDate) < datenum('01/01/2016'), :);
    bedTable16 = stationTable(and(datenum(stationTable.CollectionDate) > datenum('01/01/2016'), datenum(stationTable.CollectionDate) < datenum('01/01/2017')), :);
    bedTable18 = stationTable(datenum(stationTable.CollectionDate) > datenum('01/01/2018'), :);

    waterSampleZs =  [0.05, 0.15, 0.25, 0.5, 0.9];
    sampleZs =  [0, 0.05, 0.15, 0.25, 0.5, 0.9];

    [waterTable15ConcNWWide, waterTable15WideJitter] = categorical2jitterMat(waterTable15.concNW, waterTable15.sampleDepth, 'missing_groups', {'0.85', '0.1'}, 'positions', waterSampleZs);
    [waterTable16ConcNWWide, waterTable16WideJitter] = categorical2jitterMat(waterTable16.concNW, waterTable16.sampleDepth, 'positions', waterSampleZs);
    [waterTable18ConcNWWide, waterTable18WideJitter] = categorical2jitterMat(waterTable18.concNW, waterTable18.sampleDepth, 'positions', waterSampleZs);

    [samplesTable15d5NWWide, samplesTable15WideJitter] = categorical2jitterMat([waterTable15.gsSummNWnorm.d5; bedTable15.gsSummBedD05], ...
                                                                               [waterTable15.sampleDepth; ones(size(bedTable15,1),1)], ...
                                                                               'missing_groups', {'0.85', '0.1'}, 'positions', sampleZs);
    [samplesTable16d5NWWide, samplesTable16WideJitter] = categorical2jitterMat([waterTable16.gsSummNWnorm.d5; bedTable16.gsSummBedD05], ...
                                                                               [waterTable16.sampleDepth; ones(size(bedTable16,1),1)], ...
                                                                               'positions', sampleZs);
    [samplesTable18d5NWWide, samplesTable18WideJitter] = categorical2jitterMat([waterTable18.gsSummNWnorm.d5; bedTable18.gsSummBedD05], ...
                                                                               [waterTable18.sampleDepth; ones(size(bedTable18,1),1)], ...
                                                                               'positions', sampleZs);
    [samplesTable15d50NWWide, samplesTable15WideJitter] = categorical2jitterMat([waterTable15.gsSummNWnorm.d50; bedTable15.gsSummBedD50], ...
                                                                               [waterTable15.sampleDepth; ones(size(bedTable15,1),1)], ...
                                                                               'missing_groups', {'0.85', '0.1'}, 'positions', sampleZs);
    [samplesTable16d50NWWide, samplesTable16WideJitter] = categorical2jitterMat([waterTable16.gsSummNWnorm.d50; bedTable16.gsSummBedD50], ...
                                                                               [waterTable16.sampleDepth; ones(size(bedTable16,1),1)], ...
                                                                               'positions', sampleZs);
    [samplesTable18d50NWWide, samplesTable18WideJitter] = categorical2jitterMat([waterTable18.gsSummNWnorm.d50; bedTable18.gsSummBedD50], ...
                                                                               [waterTable18.sampleDepth; ones(size(bedTable18,1),1)], ...
                                                                               'positions', sampleZs);
    [samplesTable15d90NWWide, samplesTable15WideJitter] = categorical2jitterMat([waterTable15.gsSummNWnorm.d90; bedTable15.gsSummBedD90], ...
                                                                               [waterTable15.sampleDepth; ones(size(bedTable15,1),1)], ...
                                                                               'missing_groups', {'0.85', '0.1'}, 'positions', sampleZs);
    [samplesTable16d90NWWide, samplesTable16WideJitter] = categorical2jitterMat([waterTable16.gsSummNWnorm.d90; bedTable16.gsSummBedD90], ...
                                                                               [waterTable16.sampleDepth; ones(size(bedTable16,1),1)], ...
                                                                               'positions', sampleZs);
    [samplesTable18d90NWWide, samplesTable18WideJitter] = categorical2jitterMat([waterTable18.gsSummNWnorm.d90; bedTable18.gsSummBedD90], ...
                                                                               [waterTable18.sampleDepth; ones(size(bedTable18,1),1)], ...
                                                                               'positions', sampleZs);

    [waterTable15WLpercentWide, waterTable15WideJitter] = categorical2jitterMat(waterTable15.gsWLpercent, waterTable15.sampleDepth, 'missing_groups', {'0.85', '0.1'}, 'positions', waterSampleZs);
    [waterTable16WLpercentWide, waterTable16WideJitter] = categorical2jitterMat(waterTable16.gsWLpercent, waterTable16.sampleDepth, 'positions', waterSampleZs);
    [waterTable18WLpercentWide, waterTable18WideJitter] = categorical2jitterMat(waterTable18.gsWLpercent, waterTable18.sampleDepth, 'positions', waterSampleZs);


    jt.waterTable15ConcNWWide = waterTable15ConcNWWide; jt.waterTable15WideJitter = waterTable15WideJitter;
    jt.waterTable16ConcNWWide = waterTable16ConcNWWide; jt.waterTable16WideJitter = waterTable16WideJitter;
    jt.waterTable18ConcNWWide = waterTable18ConcNWWide; jt.waterTable18WideJitter = waterTable18WideJitter;
    jt.samplesTable15d5NWWide = samplesTable15d5NWWide; jt.samplesTable15WideJitter = samplesTable15WideJitter;
    jt.samplesTable16d5NWWide = samplesTable16d5NWWide; jt.samplesTable16WideJitter = samplesTable16WideJitter;
    jt.samplesTable18d5NWWide = samplesTable18d5NWWide; jt.samplesTable18WideJitter = samplesTable18WideJitter;
    jt.samplesTable15d50NWWide = samplesTable15d50NWWide; jt.samplesTable15WideJitter = samplesTable15WideJitter;
    jt.samplesTable16d50NWWide = samplesTable16d50NWWide; jt.samplesTable16WideJitter = samplesTable16WideJitter;
    jt.samplesTable18d50NWWide = samplesTable18d50NWWide; jt.samplesTable18WideJitter = samplesTable18WideJitter;
    jt.samplesTable15d90NWWide = samplesTable15d90NWWide; jt.samplesTable15WideJitter = samplesTable15WideJitter;
    jt.samplesTable16d90NWWide = samplesTable16d90NWWide; jt.samplesTable16WideJitter = samplesTable16WideJitter;
    jt.samplesTable18d90NWWide = samplesTable18d90NWWide; jt.samplesTable18WideJitter = samplesTable18WideJitter;
    jt.waterTable15WLpercentWide = waterTable15WLpercentWide; jt.waterTable15WideJitter = waterTable15WideJitter;
    jt.waterTable16WLpercentWide = waterTable16WLpercentWide; jt.waterTable16WideJitter = waterTable16WideJitter;
    jt.waterTable18WLpercentWide = waterTable18WLpercentWide; jt.waterTable18WideJitter = waterTable18WideJitter;

    if false
        % summary data figure
        m_plot = 0.075; % margins
        s_plot = 0.075; % spacing
        nx = 3;
        ny = 3;
        w_plot = (1-2*m_plot-(nx-1)*s_plot) / nx;
        h_plot = (1-2*m_plot-(ny-1)*s_plot) / ny;
        x_pos = repmat(m_plot+(0:(nx-1))*w_plot+(0:(nx-1))*s_plot, ny, 1)';
        y_pos = repmat(1-(h_plot+m_plot+(0:(ny-1))*h_plot+(0:(ny-1))*s_plot), nx, 1);
        exNums = [27 14 5]; % df idx to plot for 2015, 2016, 2018 data
        figure()
        for p = 1:length(exNums)
            if datenum(df(exNums(p)).CollectionDate) < datenum('01/01/2016') 
                yrCol = col1;
                yrMark = mar1;
            elseif datenum(df(exNums(p)).CollectionDate) < datenum('01/01/2017')
                yrCol = col2;
                yrMark = mar2;
            else
                yrCol = col3;
                yrMark = mar3;
            end      
            % subplot(2*length(exNums), length(exNums), [p p+length(exNums)]); hold on;
            subplot('position', [x_pos(p) y_pos(p) w_plot, h_plot]); hold on;
                % gs class specific models in background
                plot(df(exNums(p)).concProf.gsClassModel.Cs, df(exNums(p)).concProf.gsClassModel.Zs, '-', 'Color', [0.5 0.5 0.5])
                plot(concatenate_cell_vectors(df(exNums(p)).waterSamplesTable.concNWbyClass), repelem(df(exNums(p)).waterSamplesTable.sampleZ, length(df(exNums(p)).waterSamplesTable.gsClass{1})), ...
                    'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', 3, 'Marker', yrMark, 'LineStyle', 'none')
                gsclassModel_lines = plot(NaN, NaN, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', ...
                    'MarkerSize', 3, 'Marker', yrMark, 'LineStyle', '-', 'Color', [0.5 0.5 0.5]);

                % main samples and models
                NW_samples = plot(df(exNums(p)).waterSamplesTable.concNW, df(exNums(p)).waterSamplesTable.sampleZ, ...
                    'k', 'MarkerFaceColor', yrCol, 'MarkerSize', 7, 'Marker', yrMark, 'LineStyle', 'none'); % NW samples
                sumGsClass_fit = plot(df(exNums(p)).concProf.gsClassModel.CsSum, df(exNums(p)).concProf.gsClassModel.ZsSum, ...
                    'k--', 'LineWidth', 1.5); % summed class fits
                total_fit = plot(df(exNums(p)).concProf.totalModel.Cs, df(exNums(p)).concProf.totalModel.Zs, ...
                    'k:', 'LineWidth', 1.5); % single fit to total
                sumGsClass_mdl = plot(df(exNums(p)).concProf.gsClassPred.CsSum, df(exNums(p)).concProf.gsClassPred.ZsSum, ...
                    'k-', 'LineWidth', 1.5); % sum of predictions
                MYsum_mdl = plot(df(exNums(p)).concProf.MYfullPred.CsSum, df(exNums(p)).concProf.MYfullPred.Zs, ...
                    'k-.', 'LineWidth', 1.5); % wright parker prediction

                % full measured samples
                wW_samples = plot(df(exNums(p)).waterSamplesTable.conc, df(exNums(p)).waterSamplesTable.sampleZ, ...
                    'k', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 5, 'Marker', yrMark, 'LineStyle', 'none');

                box on; set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
                xlabel('conc (g/L)')
                ylabel('dist. above bed (m)')
%                 legend([NW_samples, wW_samples, total_fit, gsclassModel_lines, sumGsClass_fit, sumGsClass_mdl, MYsum_mdl], ...
%                     {'NW samples', 'full samples', 'total best-fit', 'gs class best-fits', 'gs class best-fits summed', 'gs class prediction', 'MY simulation'})
            %subplot(2*length(exNums), length(exNums), [p+length(exNums)+length(exNums) p+length(exNums)+length(exNums)+length(exNums)]);
            subplot('position', [x_pos(p+nx) y_pos(p+nx) w_plot, h_plot]); hold on;
                semilogx(str2double(df(exNums(p)).bedData.gsDistBed.Properties.RowNames), cumsum(table2array(df(exNums(p)).bedData.gsDistBed), 'omitnan'), ...
                    'k-', 'LineWidth', 2)
                hold on;
                plot(str2double(df(exNums(p)).nearBedData.gsDistNearBedNWnorm.Properties.RowNames), cumsum(table2array(df(exNums(p)).nearBedData.gsDistNearBedNWnorm), 'omitnan'), ...
                    '-', 'Color', yrCol, 'LineWidth', 2)
                xlim([1 250])
                ylim([0 100])
                box on; set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
                ylabel('% finer')
                xlabel('grain size (um)')
                legend('bed', 'ex. near-bed', 'Location', 'NorthWest')
    %         subplot(2*length(exNums), length(exNums), [p+length(exNums)+length(exNums)+length(exNums) p+length(exNums)+length(exNums)+length(exNums)+length(exNums)])
        end
        set(gcf, 'Pos', [100 100 700 700]);
        set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
        drawnow; pause(0.1);
        if printOut
            print('-dpng', '-r300', './figsExport/threeYearsExample.png');
            print('-depsc', '-painters', '-r300', './figsExport/threeYearsExample.eps');
            pause(1);
        end

        % concentration plots to fill out
        figure()
        subplot(6, 3, [1 4]); hold on;
            plot(fliplr(waterTable15ConcNWWide), waterTable15WideJitter, 'LineStyle', 'none', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
            boxplot(fliplr(waterTable15ConcNWWide), 'positions', waterSampleZs, 'orientation', 'horizontal')
                xlim([0 50])
                ylim([-0.1 1])
                ylabel('height above bed (z/H)')
                xlabel('concentration (g/L)')
                set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.50','0.9'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        subplot(6, 3, [2 5]); hold on;
            plot(fliplr(waterTable16ConcNWWide), waterTable16WideJitter, 'LineStyle', 'none', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
            boxplot(fliplr(waterTable16ConcNWWide), 'positions', waterSampleZs, 'orientation', 'horizontal')
                xlim([0 50])
                ylim([-0.1 1])
                ylabel('height above bed (z/H)')
                xlabel('concentration (g/L)')
                set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.50','0.9'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        subplot(6, 3, [3 6]); hold on;
            plot(fliplr(waterTable18ConcNWWide), waterTable18WideJitter, 'LineStyle', 'none', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
            boxplot(fliplr(waterTable18ConcNWWide), 'positions', waterSampleZs, 'orientation', 'horizontal')
                xlim([0 50])
                ylim([-0.1 1])
                ylabel('height above bed (z/H)')
                xlabel('concentration (g/L)')
                set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.50','0.9'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        subplot(6, 3, [7 10]); hold on;
            plot(fliplr(samplesTable15d5NWWide), samplesTable15WideJitter, 'LineStyle', 'none', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
            boxplot(fliplr(samplesTable15d5NWWide), 'positions', sampleZs, 'orientation', 'horizontal')
            plot(fliplr(samplesTable15d50NWWide),samplesTable15WideJitter, 'LineStyle', 'none', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
            boxplot(fliplr(samplesTable15d50NWWide), 'positions', sampleZs, 'orientation', 'horizontal')
            plot(fliplr(samplesTable15d90NWWide), samplesTable15WideJitter, 'LineStyle', 'none', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
            boxplot(fliplr(samplesTable15d90NWWide), 'positions', sampleZs, 'orientation', 'horizontal')
                xlim([0 200])
                ylim([-0.1 1])
                ylabel('height above bed (z/H)')
                xlabel('$D_{5}$, $D_{50}$, $D_{90}$ grain size $(\mu m)$')
                set(gca, 'YTick', sampleZs, 'YTickLabel', {'bed', '0.05','0.15','0.25','0.50','0.9'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        subplot(6, 3, [8 11]); hold on;
            plot(fliplr(samplesTable16d5NWWide), samplesTable16WideJitter, 'LineStyle', 'none', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
            boxplot(fliplr(samplesTable16d5NWWide), 'positions', sampleZs, 'orientation', 'horizontal')
            plot(fliplr(samplesTable16d50NWWide),samplesTable16WideJitter, 'LineStyle', 'none', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
            boxplot(fliplr(samplesTable16d50NWWide), 'positions', sampleZs, 'orientation', 'horizontal')
            plot(fliplr(samplesTable16d90NWWide), samplesTable16WideJitter, 'LineStyle', 'none', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
            boxplot(fliplr(samplesTable16d90NWWide), 'positions', sampleZs, 'orientation', 'horizontal')
                xlim([0 200])
                ylim([-0.1 1])
                ylabel('height above bed (z/H)')
                xlabel('$D_{5}$, $D_{50}$, $D_{90}$ grain size $(\mu m)$')
                set(gca, 'YTick', sampleZs, 'YTickLabel', {'bed', '0.05','0.15','0.25','0.50','0.9'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        subplot(6, 3, [9 12]); hold on;
            plot(fliplr(samplesTable18d5NWWide), samplesTable18WideJitter, 'LineStyle', 'none', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
            boxplot(fliplr(samplesTable18d5NWWide), 'positions', sampleZs, 'orientation', 'horizontal')
            plot(fliplr(samplesTable18d50NWWide),samplesTable18WideJitter, 'LineStyle', 'none', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
            boxplot(fliplr(samplesTable18d50NWWide), 'positions', sampleZs, 'orientation', 'horizontal')
            plot(fliplr(samplesTable18d90NWWide), samplesTable18WideJitter, 'LineStyle', 'none', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
            boxplot(fliplr(samplesTable18d90NWWide), 'positions', sampleZs, 'orientation', 'horizontal')
                xlim([0 200])
                ylim([-0.1 1])
                ylabel('height above bed (z/H)')
                xlabel('$D_{5}$, $D_{50}$, $D_{90}$ grain size $(\mu m)$')
                set(gca, 'YTick', sampleZs, 'YTickLabel', {'bed', '0.05','0.15','0.25','0.50','0.9'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')        
        set(gcf, 'Pos', [100 100 700 700]);
        set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
        drawnow; pause(0.1);
        if printOut
            print('-dpng', '-r300', './figsExport/threeYearsConcAndGrainSize.png');
            print('-depsc', '-painters', '-r300', './figsExport/threeYearsConcAndGrainSize.eps');
            pause(1);
        end

    %     [stationTableWLcutoffWide, stationTableWLcutoffWideJitter] = categorical2jitterMat(stationTable.gsWLcutoff, ...
    %         datenum(stationTable.CollectionDate) > datenum('01/01/2016'));
        [stationTableWLcutoffWide, stationTableWLcutoffWideJitter] = categorical2jitterMat(stationTable.gsWLcutoff, ...
            arrayfun(@(x) year(x.CollectionDate), df));
        figure(); hold on;
            plot(stationTableWLcutoffWideJitter(:,1), stationTableWLcutoffWide(:,1), 'LineStyle', 'none', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
            plot(stationTableWLcutoffWideJitter(:,2), stationTableWLcutoffWide(:,2), 'LineStyle', 'none', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
            plot(stationTableWLcutoffWideJitter(:,3), stationTableWLcutoffWide(:,3), 'LineStyle', 'none', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
            boxplot(stationTableWLcutoffWide)
    %             ylim([0 35])
                xlabel('year')
                ylabel('washload cutoff')
                set(gca, 'XTick', 1:5, 'XTickLabel', {'2015', '2016', '2018'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [100 100 350 350]);
        set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
        drawnow; pause(0.1);
        if printOut
            print('-dpng', '-r300', './figsExport/threeYearsWLcutoff.png');
            print('-depsc', '-painters', '-r300', './figsExport/threeYearsWLcutoff.eps');
            pause(1);
        end

        % make jittermat for the D50 data
        [stationTableBedD05Wide, stationTableBedD05WideJitter] = categorical2jitterMat(stationTable.gsSummBedD05, ...
            arrayfun(@(x) year(x.CollectionDate), df));
        [stationTableBedD50Wide, stationTableBedD50WideJitter] = categorical2jitterMat(stationTable.gsSummBedD50, ...
            arrayfun(@(x) year(x.CollectionDate), df));
        [stationTableBedD90Wide, stationTableBedD90WideJitter] = categorical2jitterMat(stationTable.gsSummBedD90, ...
            arrayfun(@(x) year(x.CollectionDate), df));
        stationTableGrainSizeWide = [stationTableBedD05Wide(:,1), stationTableBedD50Wide(:,1), stationTableBedD90Wide(:,1), ...
                                     stationTableBedD05Wide(:,2), stationTableBedD50Wide(:,2), stationTableBedD90Wide(:,2), ...
                                     stationTableBedD05Wide(:,3), stationTableBedD50Wide(:,3), stationTableBedD90Wide(:,3)];
        stationTableGrainSizeWideJitter = [stationTableWLcutoffWideJitter(:,1),   stationTableBedD50WideJitter(:,1)+1, stationTableBedD90WideJitter(:,1)+2, ...
                                           stationTableWLcutoffWideJitter(:,2)+2, stationTableBedD50WideJitter(:,2)+3, stationTableBedD90WideJitter(:,2)+4, ...
                                           stationTableWLcutoffWideJitter(:,3)+4, stationTableBedD50WideJitter(:,3)+5, stationTableBedD90WideJitter(:,3)+6];
        figure(); hold on;
            plot(stationTableGrainSizeWideJitter(:,(1:3)), stationTableGrainSizeWide(:,(1:3)), 'LineStyle', 'none', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
            plot(stationTableGrainSizeWideJitter(:,(4:6)), stationTableGrainSizeWide(:,(4:6)), 'LineStyle', 'none', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
            plot(stationTableGrainSizeWideJitter(:,(7:9)), stationTableGrainSizeWide(:,(7:9)), 'LineStyle', 'none', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
            boxplot(stationTableGrainSizeWide)
                ylabel('bed grain size (um)')
                set(gca, 'XTick', 1:9, 'XTickLabel', {'d5', 'd50', 'd90'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                 box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [100 100 350 350]);
        set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
        drawnow; pause(0.1);
        if printOut
            print('-dpng', '-r300', './figsExport/threeYearsBedGrainSize.png');
            print('-depsc', '-painters', '-r300', './figsExport/threeYearsBedGrainSize.eps');
            pause(1);
        end

        figure()
        subplot(1, 3, 1); hold on;
            plot(fliplr(waterTable15WLpercentWide), waterTable15WideJitter, 'o', 'Color', col1, 'MarkerSize', 4, 'Marker', markers{1})
            boxplot(fliplr(waterTable15WLpercentWide), 'positions', waterSampleZs, 'orientation', 'horizontal')
                xlim([0 100])
                ylabel('height above bed (z/H)')
                xlabel('washload percentage (%)')
                set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.5','0.9'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        subplot(1, 3, 2); hold on;
            plot(fliplr(waterTable16WLpercentWide), waterTable16WideJitter, 'o', 'Color', col2, 'MarkerSize', 4, 'Marker', markers{2})
            boxplot(fliplr(waterTable16WLpercentWide), 'positions', waterSampleZs, 'orientation', 'horizontal')
                xlim([0 100])
                ylabel('height above bed (z/H)')
                xlabel('washload percentage (%)')
                set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.5','0.9'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        subplot(1, 3, 3); hold on;
            plot(fliplr(waterTable18WLpercentWide), waterTable18WideJitter, 'o', 'Color', col3, 'MarkerSize', 4, 'Marker', markers{3})
            boxplot(fliplr(waterTable18WLpercentWide), 'positions', waterSampleZs, 'orientation', 'horizontal')
                xlim([0 100])
                ylabel('height above bed (z/H)')
                xlabel('washload percentage (%)')
                set(gca, 'YTick', waterSampleZs, 'YTickLabel', {'0.05','0.15','0.25','0.5','0.9'})
                boxes = findobj(gca, 'tag', 'Box');
                medians = findobj(gca, 'tag', 'Median');
                outliers = findobj(gca, 'tag', 'Outliers');
                set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
                set(medians, 'LineWidth', 1, 'Color', [0 0 0])
                set(outliers, 'Marker', 'none')
                box on
                set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')
        set(gcf, 'Pos', [100 100 1050 350]);
        set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off')
        drawnow; pause(0.1);
        if printOut
            print('-dpng', '-r300', './figsExport/threeYearsWLpercent.png');
            print('-depsc', '-painters', '-r300', './figsExport/threeYearsWLpercent.eps');
            pause(1);
        end
    end
    
    %% washload cutoff and perc vs ustar
    
    

    
    %% save data for export
    disp('saving.....')
    save('./dataExport/toSelectedPlots.mat', 'df', 'idx', 'nearBed', 'bed', 'station', 'mods', 'jt', 'reg', 'sdif', 'entr')
    
end
