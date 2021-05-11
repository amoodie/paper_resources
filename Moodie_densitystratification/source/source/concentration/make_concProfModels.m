function [stationStruct] = make_concProfModels(stationStruct)


    %% fit model to field data NWsamples, entire distribution at once
    stnTotalModel.params.flowDepth = stationStruct.FlowDepthAtCollection;
    stnTotalModel.params.b = stnTotalModel.params.flowDepth * 0.05;
    stnTotalModel.params.Hbb = (stnTotalModel.params.flowDepth - stnTotalModel.params.b) / stnTotalModel.params.b; % this is ((H-b)/b), a constant
    [stnTotalModel.params] = find_RouseModel(stnTotalModel.params, stationStruct.waterSamplesTable.sampleDepth, ...
        stationStruct.waterSamplesTable.sampleZ, stationStruct.waterSamplesTable.concNW);
    [stnTotalModel.Zs, stnTotalModel.Cs] = evaluate_RouseModel(stnTotalModel.params.flowDepth, stnTotalModel.params.cb, stnTotalModel.params.Hbb, stnTotalModel.params.Rou);
    [stnTotalModel.ZsNorm, stnTotalModel.CsNorm] = normalize_model(stnTotalModel.Zs, stnTotalModel.Cs);
    stnTotalModel.Cbar = trapz(stnTotalModel.Zs, stnTotalModel.Cs) / stnTotalModel.params.flowDepth;
    stationStruct.concProf.totalModel = stnTotalModel;
    
    
    %% fit model to field data gs_specific models, and then combine proportionately
    for c = 1:size(stationStruct.nearBedData.gsDistNearBedNWnorm,1)
        stnGsClassModelc(c).gsClass = mean(arrayfun(@(x) x{:}(c), stationStruct.waterSamplesTable.gsClass)); % stationStruct.waterSamplesTable.gsClass{1}(c);
        stnGsClassModelc(c).params.flowDepth = stationStruct.FlowDepthAtCollection;
        stnGsClassModelc(c).params.b = stnGsClassModelc(c).params.flowDepth * 0.05;
        stnGsClassModelc(c).params.Hbb = (stnGsClassModelc(c).params.flowDepth - stnGsClassModelc(c).params.b) / stnGsClassModelc(c).params.b; % this is ((H-b)/b), a constant
        [stnGsClassModelc(c).params] = find_RouseModel(stnGsClassModelc(c).params, stationStruct.waterSamplesTable.sampleDepth, ...
            stationStruct.waterSamplesTable.sampleZ, ...
            arrayfun(@(x) x{:}(c), stationStruct.waterSamplesTable.concNWbyClass));
        paramsStruct(c) = stnGsClassModelc(c).params;
        [stnGsClassModelc(c).Zs, stnGsClassModelc(c).Cs] = evaluate_RouseModel(stnGsClassModelc(c).params.flowDepth, stnGsClassModelc(c).params.cb, stnGsClassModelc(c).params.Hbb, stnGsClassModelc(c).params.Rou);
        [stnGsClassModelc(c).ZsNorm, stnGsClassModelc(c).CsNorm] = normalize_model(stnGsClassModelc(c).Zs, stnGsClassModelc(c).Cs);
        stnGsClassModelc(c).Cbar = trapz(stnGsClassModelc(c).Zs, stnGsClassModelc(c).Cs) / stnGsClassModelc(c).params.flowDepth;
        stationStruct.concProf.gsClassModel.gsClass(c) = stnGsClassModelc(c);
    end
    stnGsClassModel.gsClass = horzcat(stnGsClassModelc.gsClass);
    stnGsClassModel.Zs = horzcat(stnGsClassModelc.Zs);
    stnGsClassModel.Cs = horzcat(stnGsClassModelc.Cs);
    stnGsClassModel.ZsSum = nanmean(stnGsClassModel.Zs, 2);
    stnGsClassModel.CsSum = nansum(stnGsClassModel.Cs, 2);
    [stnGsClassModel.ZsNorm, stnGsClassModel.CsNorm] = normalize_model(stnGsClassModel.ZsSum, stnGsClassModel.CsSum);
    stnGsClassModel.params = paramsStruct;
    stnGsClassModel.Rou = horzcat(stnGsClassModel.params.Rou);
    % stnGsClassModel.bulkRou = nansum(vertcat(stnGsClassModel.params.Rou) .* (table2array(stationStruct.nearBedData.gsDistNearBedNWnorm) ./ 100)); % by weighted avg
    stnGsClassModel.RouFix = stnGsClassModel.Rou;
    stnGsClassModel.bulkRou = trapz(cumsum(table2array(stationStruct.nearBedData.gsDistNearBedNWnorm) ./ 100)', stnGsClassModel.RouFix); % bulk by integration
    stationStruct.concProf.gsClassModel = stnGsClassModel;
    
    %% make prediction model from the data
    con = load_conset('quartz-water');
    stnPred.params.flowDepth = stationStruct.FlowDepthAtCollection;
    stnPred.params.ustar = stationStruct.Velocity.ustarCalib;
    stnPred.params.b = stnPred.params.flowDepth * 0.05;
    stnPred.params.Hbb = (stnPred.params.flowDepth - stnPred.params.b) / stnPred.params.b; % this is ((H-b)/b), a constant
    stnPred.params.cb = mean(stationStruct.waterSamplesTable.concNW(stationStruct.waterSamplesTable.sampleDepth == 0.95)) ...
        .* (table2array(stationStruct.nearBedData.gsDistNearBedNWnorm) ./ 100)';
    stnPred.nearBedClass = cellfun(@(x) str2double(x), stationStruct.nearBedData.gsDistNearBedNWnorm.Properties.RowNames)';
    stnPred.nearBedClassMid = stnPred.nearBedClass; % [NaN ((stnPred.nearBedClass(2:end) -  stnPred.nearBedClass(1:end-1)) / 2) +  stnPred.nearBedClass(1:end-1)];
    stnPred.nearBedClassPercs = (table2array(stationStruct.nearBedData.gsDistNearBedNWnorm) ./ 100)';
    stnPred.params.ws = get_DSV(stnPred.nearBedClassMid .* 1e-6, 0.7, 3.5, con);
    stnPred.params.Rou = stnPred.params.ws ./ (0.41 .* stnPred.params.ustar);
    for c = 1:length(stnPred.nearBedClassMid)
        [stnPred.Zs(:, c), stnPred.Cs(:, c)] = evaluate_RouseModel(stnPred.params.flowDepth, ...
            stnPred.params.cb(c), stnPred.params.Hbb, stnPred.params.Rou(c));
        [stnPred.ZsNorm(:, c), stnPred.CsNorm(:, c)] = normalize_model(stnPred.Zs(:, c), stnPred.Cs(:, c));
    end
    stnPred.CsSum = nansum(stnPred.Cs, 2);
    stnPred.ZsSum = nanmean(stnPred.Zs, 2);
    [stnPred.ZsSumNorm, stnPred.CsSumNorm] = normalize_model(stnPred.ZsSum, stnPred.CsSum);
    stnPred.Rou = stnPred.params.Rou;
    stnPred.bulkRou = nansum(stnPred.params.Rou .* (table2array(stationStruct.nearBedData.gsDistNearBedNWnorm) ./ 100)'); % bulk by weighted avg
    % stnPred.bulkRou = trapz(cumsum(table2array(stationStruct.nearBedData.gsDistNearBedNWnorm) ./ 100)', stnPred.params.Rou); % bulk by integration
    stnPred.d50Rou = get_DSV(stationStruct.nearBedData.gsSummNearBedNWnorm.d50*1e-6, 0.7, 3.5, con) / ...
        (0.41 * stnPred.params.ustar);
    stnPred.Cbari = trapz(stnPred.Zs(:,1), stnPred.Cs, 1) ./ stnTotalModel.params.flowDepth;
    stnPred.CbariSum = nansum(stnPred.Cbari, 2);
    stationStruct.concProf.gsClassPred = stnPred;
    
    %% make a profile using the alpha given by WP04
    stnWPalpha.params = stnPred.params;
    stnWPalpha.params.alpha = calculate_alpha_WP04(nansum(stnWPalpha.params.cb)./2650, stationStruct.Velocity.slope);
    stnWPalpha.params.Rou = stnWPalpha.params.Rou .* (1./stnWPalpha.params.alpha);
    for c = 1:length(stnWPalpha.params.Rou)
        [stnWPalpha.Zs(:, c), stnWPalpha.Cs(:, c)] = evaluate_RouseModel(stnWPalpha.params.flowDepth, ...
            stnWPalpha.params.cb(c), stnWPalpha.params.Hbb, stnWPalpha.params.Rou(c));
        [stnWPalpha.ZsNorm(:, c), stnWPalpha.CsNorm(:, c)] = normalize_model(stnWPalpha.Zs(:, c), stnWPalpha.Cs(:, c));
    end
    stnWPalpha.CsSum = nansum(stnWPalpha.Cs, 2);
    stnWPalpha.ZsSum = nanmean(stnWPalpha.Zs, 2);
    [stnWPalpha.ZsSumNorm, stnWPalpha.CsSumNorm] = normalize_model(stnWPalpha.ZsSum, stnWPalpha.CsSum);
    stnWPalpha.alpha = stnWPalpha.params.alpha;
    stationStruct.concProf.WPalphaPred = stnWPalpha;
    
    %% make a plot
    if false
        figure(); hold on;
        set(gca,'ColorOrderIndex',1)
        %% dashes for the fits
        plot(stationStruct.concProf.gsClassModel.Cs, stationStruct.concProf.gsClassModel.Zs, '--')
        plot(stationStruct.concProf.gsClassModel.CsSum, stationStruct.concProf.gsClassModel.ZsSum, 'k--')
        set(gca,'ColorOrderIndex',1)
        %% solids for the preds
        plot(stationStruct.concProf.gsClassPred.Cs, stationStruct.concProf.gsClassPred.Zs, '-')
        plot(stationStruct.concProf.gsClassPred.CsSum, stationStruct.concProf.gsClassPred.ZsSum, 'k-')
    end
    
end