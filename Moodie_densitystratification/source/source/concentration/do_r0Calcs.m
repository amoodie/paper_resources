function [r0calcs] = do_r0Calcs(stationStruct, modelFlag)

    if strcmp(modelFlag, 'total')
        concProfModel = stationStruct.concProf.totalModel;
        concProfModel.CsSum = concProfModel.Cs;
        concProfModel.ZsSum = concProfModel.Zs;
    elseif strcmp(modelFlag, 'gsClass')
        concProfModel = stationStruct.concProf.gsClassModel;
    end

    %% handler for the various velocity profiles coming back from processing
    % if velocity profile exists, get list, otherwise placeholder for skip
    % to 'else' below
    if ~isempty(stationStruct.velProf)
        velProfs_list = fieldnames(stationStruct.velProf);
    else
        velProfs_list = {''};
    end
    
    % check for each possible comeback
    if ismember('modelAlpha', velProfs_list)
        r0calcs.Uprof = interp1(stationStruct.velProf.modelAlpha.Zs, ...
            stationStruct.velProf.modelAlpha.Us, concProfModel.ZsSum);
    elseif ismember('modelNoAlpha', velProfs_list)
        r0calcs.Uprof = interp1(stationStruct.velProf.modelNoAlpha.Zs, ...
            stationStruct.velProf.modelNoAlpha.Us, concProfModel.ZsSum);
    else
        r0calcs.Uprof = (stationStruct.Velocity.ustarCalib / (0.41)) .* log(concProfModel.ZsSum ./ (3 * stationStruct.bedData.gsSummBed.d90 * 1e-6 / 30));
    end
    
    %% do calcualtions
    r0calcs.udz = trapz(concProfModel.ZsSum, r0calcs.Uprof);
    r0calcs.Ubar = r0calcs.udz ./ stationStruct.FlowDepthAtCollection; % depth avg vel
    r0calcs.ucdz = trapz(concProfModel.ZsSum, (r0calcs.Uprof .* concProfModel.CsSum));
    r0calcs.Cm = (r0calcs.ucdz ./ r0calcs.udz) ./ 2650;
    r0calcs.Cbar = r0calcs.ucdz / (stationStruct.FlowDepthAtCollection * r0calcs.Ubar); % depth flux averaged concentration
    r0calcs.cb = concProfModel.CsSum(1);
    r0calcs.r0 = r0calcs.cb ./ r0calcs.Cbar;
     
end