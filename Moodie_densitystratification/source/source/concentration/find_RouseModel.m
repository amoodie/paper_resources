function [modelParams] = find_RouseModel(modelParams, sampleDepth, sampleZ, conc)
    % receives modelParams.b only!
    samplesTable = array2table([sampleDepth, sampleZ, conc], 'VariableNames', {'sampleDepth', 'sampleZ', 'conc'});
    
    % using try catch statements below because sometimes it fails due to bad data,
    % and sometimes I want to force an error for bad values
    
    % set some boundaries for the fit with more free parameters, to ensure good convergence in the end
    RouMin = 0;
    RouMax = 5;
    cbMin = nanmean(samplesTable.conc(samplesTable.sampleDepth > 0.9)) - (2*std(samplesTable.conc(samplesTable.sampleDepth > 0.9), 'omitnan'));
    cbMax = nanmean(samplesTable.conc(samplesTable.sampleDepth > 0.9)) + (2*std(samplesTable.conc(samplesTable.sampleDepth > 0.9), 'omitnan'));
    
    % try first with cb and Rou as free parameters
    nTry = 10;
    goodFit = false;
    cnt = 0;
    RouGuess = 0.5;
    while ~or(cnt>nTry, goodFit) % while cnt>try or no good fit
        try
            modelParams.formula = strcat('conc ~ cb * ( ((', num2str(modelParams.flowDepth), '-sampleZ)/sampleZ) / ', num2str(modelParams.Hbb), ' )^Rou');
            modelParams.type = 'a';
            modelParams.guess = [nanmean(samplesTable.conc(samplesTable.sampleDepth > 0.9)) RouGuess]; % guesses for cb and Rou
            opts = statset('Robust', 'on', 'RobustWgtFun', 'bisquare');
            modelParams.model = fitnlm(samplesTable, modelParams.formula, modelParams.guess, 'Options', opts);
            modelParams.Rou = modelParams.model.Coefficients.Estimate(1);
            modelParams.cb = modelParams.model.Coefficients.Estimate(2);

            if and(and(modelParams.Rou > RouMin, modelParams.Rou < RouMax), and(modelParams.cb > cbMin, modelParams.cb < cbMax))
                goodFit = true;
            else
                cnt = cnt + 1;
                RouGuess = RouGuess .* (rand(1)+0.5);
            end
        catch
            cnt = cnt + 1;
        end
    end

    
    % if fails, fix cb and find only Rou
    if ~goodFit
    
        cnt = 0;
        while ~or(cnt>nTry, goodFit) % while cnt>try or no good fit
            try
                modelParams.cb = nanmean(samplesTable.conc(samplesTable.sampleDepth > 0.9));
                modelParams.formula = strcat('conc ~ ', num2str(modelParams.cb, '%.15f'),' * ( ((', num2str(modelParams.flowDepth), '-sampleZ)/sampleZ) / ', num2str(modelParams.Hbb), ' )^Rou');
                modelParams.type = 'b';
                modelParams.guess = 0.5; % guess for Rou
                modelParams.model = fitnlm(samplesTable, modelParams.formula, modelParams.guess);
                modelParams.Rou = modelParams.model.Coefficients.Estimate(1);

                if and(modelParams.Rou > 0, modelParams.Rou < 5)
                    goodFit = true;
                else
                    cnt = cnt + 1;
                end
            catch
                cnt = cnt + 1;
            end
        end
    end


    % if this fails, assume Rou=0, cb=mean(C(z))
    if ~goodFit        
        warning('failed to determine model fit of Type A or B; assume Rou=0, cb=mean(C(z)')
        modelParams.cb = mean(conc);
        modelParams.formula = NaN;
        modelParams.type = NaN;
        modelParams.guess = NaN; % guess for Rou
        modelParams.model = NaN;
        modelParams.Rou = 0;
    end
    
    
end