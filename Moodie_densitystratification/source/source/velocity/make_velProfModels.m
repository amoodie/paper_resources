function [stationStruct] = make_velProfModels(stationStruct)
    % main loop for velocity profile models for each station
    %   create one velocity profile model for each station, 
    %   using the combined profile model
    
    modelParams.ks = 3 * stationStruct.bedData.gsSummBed.d90 * 1e-6;
    modelParams.z0 = modelParams.ks;

    if ~isempty(stationStruct.Velocity.data.measZ)
        
        stationStruct.velProf.hasData = true;
        
        % make and evaluate a model which has no alpha parameter
        try
            [stationStruct.velProf.modelNoAlpha.params] = find_LotwModel(modelParams, stationStruct.Velocity.data.measZ, ...
                stationStruct.Velocity.data.velMean, 'fix');
        catch
            [stationStruct.velProf.modelNoAlpha.params] = find_LotwModel(modelParams, stationStruct.Velocity.data.measZ(~stationStruct.Velocity.data.measZ==0), ...
                stationStruct.Velocity.data.velMean(~stationStruct.Velocity.data.measZ==0), 'fix');
        end
        [stationStruct.velProf.modelNoAlpha.Zs, stationStruct.velProf.modelNoAlpha.Us] = evaluate_LotwModel(stationStruct.FlowDepthAtCollection, ...
            stationStruct.velProf.modelNoAlpha.params.z0, stationStruct.velProf.modelNoAlpha.params.ustar, stationStruct.velProf.modelNoAlpha.params.alpha);
        
        % make and evaluate a model which has an alpha parameter
        [stationStruct.velProf.modelAlpha.params] = find_LotwModel(modelParams, stationStruct.Velocity.data.measZ, ...
            stationStruct.Velocity.data.velMean, 'free');
        [stationStruct.velProf.modelAlpha.Zs, stationStruct.velProf.modelAlpha.Us] = evaluate_LotwModel(stationStruct.FlowDepthAtCollection, ...
            stationStruct.velProf.modelAlpha.params.z0, stationStruct.velProf.modelAlpha.params.ustar, stationStruct.velProf.modelAlpha.params.alpha);
        
        % make a model to determine the regression in log-linear space.
        [stationStruct.velProf.modelKred.params] = find_LotwModel(modelParams, stationStruct.Velocity.data.measZ, ...
            stationStruct.Velocity.data.velMean, 'log-linear', stationStruct.Velocity.ustarCalib);
        [stationStruct.velProf.modelKred.Zs, stationStruct.velProf.modelKred.Us] = evaluate_LotwModel(stationStruct.FlowDepthAtCollection, ...
            stationStruct.velProf.modelKred.params.z0, stationStruct.velProf.modelKred.params.ustar, stationStruct.velProf.modelKred.params.alpha);

        
        if false
        
            figure(); hold on;
            title(stationStruct.StationID, 'Interpreter', 'none')
            if ~isempty(stationStruct.Velocity.meter)
                plot(stationStruct.Velocity.meterTable.mean, stationStruct.Velocity.meterTable.measZ, 'ok', 'MarkerSize', 5)
                plot(stationStruct.velProf.modelMeterNoAlpha.Us, stationStruct.velProf.modelMeterNoAlpha.Zs, 'r--')
                plot(stationStruct.velProf.modelMeterAlpha.Us, stationStruct.velProf.modelMeterAlpha.Zs, 'r-')
            end
            if ~isempty(stationStruct.Velocity.adcp)
                plot(stationStruct.Velocity.adcpTable.mean, stationStruct.Velocity.adcpTable.measZ, '^k', 'MarkerSize', 5)
                plot(stationStruct.velProf.modelAdcpNoAlpha.Us, stationStruct.velProf.modelAdcpNoAlpha.Zs, 'g--')
                plot(stationStruct.velProf.modelAdcpAlpha.Us, stationStruct.velProf.modelAdcpAlpha.Zs, 'g-')
            end
        
        end
        
       
    else
        stationStruct.velProf.hasData = false;
    end
    
    
    %% predictive models
    % evaluate a no stratification log-law model given the ustar
    [stationStruct.velProf.noStratLoglaw.Zs, stationStruct.velProf.noStratLoglaw.Us] = evaluate_LotwModel(stationStruct.FlowDepthAtCollection, ...
            modelParams.z0, stationStruct.Velocity.ustarCalib, 1);
        
    % make a profile using the alpha given by WP04
    cb = mean(stationStruct.waterSamplesTable.concNW(stationStruct.waterSamplesTable.sampleDepth == 0.95)) ...
        .* (table2array(stationStruct.nearBedData.gsDistNearBedNWnorm) ./ 100)';
    stnWPalpha.alpha = calculate_alpha_WP04(nansum(cb)./2650, stationStruct.Velocity.slope);
    [stnWPalpha.Zs, stnWPalpha.Us] = evaluate_LotwModel(stationStruct.FlowDepthAtCollection, ...
            modelParams.z0, stationStruct.Velocity.ustarCalib, stnWPalpha.alpha);
    stationStruct.velProf.WPalphaPred = stnWPalpha;
end
