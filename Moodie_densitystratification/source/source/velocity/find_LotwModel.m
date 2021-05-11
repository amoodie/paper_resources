function [modelParams] = find_LotwModel(modelParams, measZ, vel, alphaSwitch, varargin)

    measTable = array2table([measZ, vel], 'VariableNames', {'measZ', 'vel'});

    switch alphaSwitch
        
        % no alpha parameter (fixed to 1)
        case 'fix'
            modelParams.alpha = 1;
            modelParams.formula = strcat('vel ~ (ustar / (', num2str(modelParams.alpha, '%f'), ' * 0.41)) * log(measZ / (', num2str(modelParams.z0, '%f'), '))');
            modelParams.model = fitnlm(measTable, modelParams.formula, 0.1);
            modelParams.ustar = modelParams.model.Coefficients.Estimate(1);
            modelParams.Rsq = modelParams.model.Rsquared.Ordinary;
            
        % allow an alpha parameter to be fit
        case 'free'
            modelParams.formula = strcat('vel ~ (ustar / (alpha * 0.41)) * log(measZ / (', num2str(modelParams.z0, '%f'), '))');
            modelParams.model = fitnlm(measTable, modelParams.formula, [0.2, 0.1]);
            modelParams.ustar = modelParams.model.Coefficients.Estimate(2);
            modelParams.alpha = modelParams.model.Coefficients.Estimate(1);
            modelParams.Rsq = modelParams.model.Rsquared.Ordinary;
            
        % tranform the y axis to log, fit linear regression
        case 'log-linear'
            measTable.logmeasZ = log(measTable.measZ / modelParams.z0);
            measTable.normVel = measTable.vel / varargin{1};
            modelParams.formula = 'normVel ~ logmeasZ -1';
            modelParams.model = fitlm(measTable, modelParams.formula);
            modelParams.kappa_alpha = modelParams.model.Coefficients.Estimate(1) ^ -1;
            modelParams.ustar = varargin{1};
            modelParams.alpha = modelParams.kappa_alpha / 0.41;
            modelParams.Rsq = modelParams.model.Rsquared.Ordinary;
            
    end
end