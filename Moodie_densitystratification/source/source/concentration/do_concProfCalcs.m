function [stationStruct] = do_concProfCalcs(stationStruct, modelFlag)
    

    %% make concentation profile calculations for each station
    % these are basically all the statistics that will be plotted in ExploreDataTable
    
    % copy the concProf structure out to operate directly on it
    i_concProf = stationStruct.concProf;
    i_velProf = stationStruct.velProf;

    %% add a few attributes to the concProf structure
    i_concProf.nearBedConcNWmodelMean = mean([i_concProf.gsClassPred.CsSum(1), i_concProf.totalModel.Cs(1)]); % best-fit average near bed concentration -- this is my BEST estimate
    i_concProf.nearBedConcNWmodelMeanGsClass = (  i_concProf.gsClassPred.Cs(1,:) + i_concProf.gsClassModel.Cs(1,:)  ) ./ 2;
    i_concProf.D50ws = get_DSV(stationStruct.nearBedData.gsSummNearBedNWnorm.d50*1e-6, 0.7, 3.5, load_conset('quartz-water')); % D50 settling velocity
    
    %% calculate some variations of the WrightParker alpha number
    i_concProf.alpha.byModelOnly = i_concProf.D50ws ./ (stationStruct.Velocity.ustarCalib .* 0.41 .* i_concProf.totalModel.params.Rou);
    
    % use cumultive prediction to compare to models
    i_concProf.alpha.byPredD50RouToTotalModelRou = i_concProf.gsClassPred.d50Rou ./ i_concProf.totalModel.params.Rou; % using d50 rouse number
    i_concProf.alpha.byPredD50RouToGsClassModelRou = i_concProf.gsClassPred.d50Rou ./ vertcat(i_concProf.gsClassModel.bulkRou);
    i_concProf.alpha.byPredBulkRouToTotalModelRou = i_concProf.gsClassPred.bulkRou ./ i_concProf.totalModel.params.Rou; % using bulk (weighted) rouse number
    i_concProf.alpha.byPredBulkRouToGsClassModelRou = i_concProf.gsClassPred.bulkRou ./ vertcat(i_concProf.gsClassModel.bulkRou);
    
    % do the calculation for each grain size class
    i_concProf.alpha.byPredGsClassRouToGsClassModelRou = i_concProf.gsClassPred.Rou ./ i_concProf.gsClassModel.Rou;
    
    i_concProf.alpha.byResid = nansum(   ( i_concProf.gsClassPred.params.ws ./ (i_concProf.totalModel.params.Rou .* 0.41 .* i_concProf.gsClassPred.params.ustar) ) .* i_concProf.gsClassPred.nearBedClassPercs   );
    
    % calculate alpha directly from the k_red of the WP model
    i_concProf.alpha.byKredMY = i_concProf.MYfullPred.K_red;
    i_concProf.alpha.byWP04func = i_concProf.WPalphaPred.alpha;
    
    %% calculate the mean signed deviation statistics
    % calculate the offset for each depth, for the summed models
    i_concProf.theta.byDepthToTotalModel = [i_concProf.gsClassPred.CsSum(1:end-1) - i_concProf.totalModel.Cs(1:end-1); NaN]; % exclude last because always zero!
    i_concProf.theta.byDepthToGsClassModel = [i_concProf.gsClassPred.CsSum(1:end-1) - i_concProf.gsClassModel.CsSum(1:end-1); NaN]; % exclude last because always zero!
    
    % calculate the offset for each depth, for the individual grain size models
    i_concProf.theta.byDepthGsClassToGsClassModel = [i_concProf.gsClassPred.Cs(1:end-1,:) - i_concProf.gsClassModel.Cs(1:end-1,:); NaN(1, size(i_concProf.gsClassPred.Cs,2))]; % exclude last because always zero!
    i_concProf.theta.byDepthMYfullGsClassToGsClassModel = [i_concProf.MYfullPred.Cs(1:end-1,:) - i_concProf.gsClassModel.Cs(1:end-1,:); NaN(1, size(i_concProf.gsClassPred.Cs,2))]; % exclude last because always zero!
    
    i_concProf.theta.byDepthGsClassToGsClassModelNorm = i_concProf.theta.byDepthGsClassToGsClassModel ./ i_concProf.gsClassModel.Cs(1,:);
    i_concProf.theta.byDepthMYfullGsClassToGsClassModelNorm = i_concProf.theta.byDepthMYfullGsClassToGsClassModel ./ i_concProf.gsClassModel.Cs(1,:);
    
    % calculate the mean deviation over the depth, for the summed models
    i_concProf.theta.meanToTotalModel = nanmean(i_concProf.theta.byDepthToTotalModel);
    i_concProf.theta.meanToGsClassModel = nanmean(i_concProf.theta.byDepthToGsClassModel);
    
    % calculate the mean deviation over the depth, for the individual grain size models
    i_concProf.theta.meanGsClassToGsClassModel = nanmean(i_concProf.theta.byDepthGsClassToGsClassModel,1);
    
    % calculate the normalized mean deviation, for the summed models
    i_concProf.theta.meanNormToTotalModel = i_concProf.theta.meanToTotalModel ./ i_concProf.nearBedConcNWmodelMean;
    i_concProf.theta.meanNormToGsClassModel = i_concProf.theta.meanToGsClassModel ./ i_concProf.nearBedConcNWmodelMean;
    
    % calculate the normalized mean deviation, for the individual grain size models
    i_concProf.theta.meanNormGsClassToGsClassModel = i_concProf.theta.meanGsClassToGsClassModel ./ i_concProf.nearBedConcNWmodelMeanGsClass;
    
    i_concProf.theta.RMSDToTotalModel = sqrt(nansum((i_concProf.theta.byDepthToTotalModel).^2, 1));
    i_concProf.theta.RMSDToGsClassModel = sqrt(nansum((i_concProf.theta.byDepthToGsClassModel).^2, 1));
    
    i_concProf.theta.stdToTotalModel = std(i_concProf.theta.byDepthToTotalModel, [], 1, 'omitnan');
    i_concProf.theta.stdToGsClassModel = std(i_concProf.theta.byDepthToGsClassModel, [], 1, 'omitnan');
    
    
    %% perform a two sample KS test to compare the distributions
    [i_concProf.kst.ToTotalModel] = kstest2(i_concProf.gsClassPred.CsSumNorm, i_concProf.totalModel.CsNorm);
    [i_concProf.kst.ToGsClassModel] = kstest2(i_concProf.gsClassPred.CsSumNorm, i_concProf.gsClassModel.CsNorm);
    
    
    %% do alpha to determine the misfit with WP model and 
    i_concProf.beta.MYr0 = i_concProf.MYfullPred.Cs(1,:) ./ mean(i_concProf.MYfullPred.Cs, 1);
    i_concProf.beta.fieldr0 = i_concProf.gsClassModel.Cs(1,:) ./ mean(i_concProf.gsClassModel.Cs, 1);
    i_concProf.beta.fieldr0Bulk = i_concProf.gsClassModel.CsSum(1,:) ./ mean(i_concProf.gsClassModel.CsSum, 1);
    i_concProf.beta.Predr0 = i_concProf.gsClassPred.Cs(1,:) ./ mean(i_concProf.gsClassPred.Cs, 1);
    
    
    %% calculate the RMSE for each of the models against the data
    [i_concProf.gsClassPred] = calculate_model_RMSE(stationStruct.waterSamplesTable, i_concProf.gsClassPred);
    [i_concProf.WPalphaPred] = calculate_model_RMSE(stationStruct.waterSamplesTable, i_concProf.WPalphaPred);
    [i_concProf.MYfullPred] = calculate_model_RMSE(stationStruct.waterSamplesTable, i_concProf.MYfullPred);
    [i_concProf.gsClassModel] = calculate_model_RMSE(stationStruct.waterSamplesTable, i_concProf.gsClassModel);
    [i_concProf.totalModel] = calculate_model_RMSE(stationStruct.waterSamplesTable, i_concProf.totalModel);
    
    
    %% sediment transport rate calculations 
    top_bins = 1:floor(size(i_concProf.gsClassModel.Cs,1)*0.2); 
    if stationStruct.velProf.hasData
        i_velProf.modelKred.UsOnCsGrid = interp1(i_velProf.modelKred.Zs, i_velProf.modelKred.Us, i_concProf.gsClassModel.ZsSum);
        i_concProf.transport.fits = trapz(i_concProf.gsClassModel.ZsSum, i_concProf.gsClassModel.CsSum .* i_velProf.modelKred.UsOnCsGrid);
        i_concProf.transport.fits_upper = trapz(i_concProf.gsClassModel.ZsSum(top_bins, :), i_concProf.gsClassModel.CsSum(top_bins, :) .* i_velProf.modelKred.UsOnCsGrid(top_bins, :));
    end
    i_velProf.noStratLoglaw.UsOnCsGrid = interp1(i_velProf.noStratLoglaw.Zs, i_velProf.noStratLoglaw.Us, i_concProf.gsClassPred.ZsSum);
    i_concProf.transport.noStrat = trapz(i_concProf.gsClassPred.ZsSum, i_concProf.gsClassPred.CsSum .* i_velProf.noStratLoglaw.UsOnCsGrid);
    i_concProf.transport.noStrat_upper = trapz(i_concProf.gsClassPred.ZsSum(top_bins, :), i_concProf.gsClassPred.CsSum(top_bins, :) .* i_velProf.noStratLoglaw.UsOnCsGrid(top_bins, :));
    i_concProf.transport.MYmodelGsClass = trapz(i_concProf.MYfullPred.ZsSum, i_concProf.MYfullPred.Cs .* i_velProf.MYfullPred.Us);
    i_concProf.transport.MYmodel = trapz(i_concProf.MYfullPred.ZsSum, i_concProf.MYfullPred.CsSum .* i_velProf.MYfullPred.Us);
    i_concProf.transport.MYmodelNoSed = trapz(i_concProf.MYfullPred.ZsSum, i_concProf.MYfullPred.CsSum_nosed .* i_velProf.MYfullPred.Us_nosed);
    ws = get_DSV(gsDistClass2num(stationStruct.nearBedData.gsDistNearBedNWnorm) ./ 1e6, 0.7, 3.5, load_conset('quartz-water'));
    i_concProf.transport.Zeff = trapz(cumsum(table2array(stationStruct.nearBedData.gsDistNearBedNWnorm) ./ 100)', ( stationStruct.Velocity.ustarCalib ./ ws ));
    i_concProf.transport.ZReff = trapz(cumsum(table2array(stationStruct.nearBedData.gsDistNearBedNWnorm) ./ 100)', ( ws ./ (0.41 .* stationStruct.Velocity.ustarCalib) ));
    i_velProf.WPalphaPred.UsOnCsGrid = interp1(i_velProf.WPalphaPred.Zs, i_velProf.WPalphaPred.Us, i_concProf.WPalphaPred.ZsSum);
    i_concProf.transport.WPalpha = trapz(i_concProf.WPalphaPred.ZsSum, i_concProf.WPalphaPred.CsSum .* i_velProf.WPalphaPred.UsOnCsGrid);
    i_concProf.transport.field_bulkRou = vertcat(i_concProf.gsClassModel.bulkRou);
    i_concProf.transport.field_Rou = i_concProf.gsClassModel.Rou;
    i_concProf.transport.pred_Rou = i_concProf.gsClassPred.Rou;
    
    % transport in the top of the water column, compute and compare D50s
    fit_concperbin = nansum(i_concProf.gsClassModel.Cs(top_bins, :),1);
    pred_concperbin = nansum(i_concProf.gsClassPred.Cs(top_bins, :),1);
    i_concProf.transport.fit_top_cumdist =  cumsum(fit_concperbin ./ nansum(fit_concperbin), 'omitnan');
    i_concProf.transport.pred_top_cumdist =  cumsum(pred_concperbin ./ nansum(pred_concperbin), 'omitnan');
    i_concProf.transport.fit_D_ii = grainSizeInterp(i_concProf.transport.fit_top_cumdist, ...
        i_concProf.gsClassModel.gsClass, [0.05, 0.50, 0.90]);
    i_concProf.transport.pred_D_ii = grainSizeInterp(i_concProf.transport.pred_top_cumdist, ...
        i_concProf.gsClassPred.nearBedClass, [0.10, 0.50, 0.90]);
    i_concProf.transport.fit_perconcbin = fit_concperbin;
    i_concProf.transport.pred_perconcbin = pred_concperbin;
    i_concProf.transport.concRatio = nansum(fit_concperbin)/nansum(pred_concperbin);
    i_concProf.transport.concRatioGsClass = fit_concperbin ./ pred_concperbin;
    i_concProf.transport.concFracReduct = (nansum(pred_concperbin)-nansum(fit_concperbin))/nansum(pred_concperbin);
    
    %% replace into original slot in stationStruct
    stationStruct.concProf = i_concProf;
    stationStruct.velProf = i_velProf;
    
end