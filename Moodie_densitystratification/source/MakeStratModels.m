function MakeStratModels()

    source_path = genpath('source');
    data_path = genpath('dataSrc');
    export_path = genpath('dataExport');
    addpath(source_path, data_path, export_path);

    %% turn off some warnings
    warning_list = {'stats:nlinfit:IllConditionedJacobian', ...
                    'MATLAB:rankDeficientMatrix', ...
                    'stats:nlinfit:ModelConstantWRTParam', ...
                    'stats:nlinfit:IterationLimitExceeded'};
    for w = 1:length(warning_list)
        warning('off', warning_list{w});
    end
    
    %% load data
    load('./dataExport/stationSurveyDataTable.mat', 'df')
    [con] = load_conset('quartz-water');
    
    %% Rouse and velocity models and calcs for each station
    for i = 1:size(df, 1)
        
        disp(['beginning station ', num2str(i), ' of ' num2str(size(df,1))])
        df(i).concProf = [];
        df(i).velProf = [];
        
        disp(['    processing velocity models....'])
        [df(i)] = make_velProfModels(df(i));
        [df(i)] = do_velProfCalcs(df(i));
        
        disp(['    processing concentration models....'])
        disp(['        regression models....'])
        [df(i)] = make_concProfModels(df(i)); % make rouse models to determine rouse number and cb=
        disp(['        WP04 model....'])
        [df(i)] = MellorYamadaConcVelModel_struct(df(i)); % adds fields to concProf and VelProf fields
        disp(['        profile calculations....'])
        modelFlag = 'gsClass';
        [df(i)] = do_concProfCalcs(df(i), modelFlag);
        [df(i).r0calcs] = do_r0Calcs(df(i), modelFlag);
        disp(['completed station ', num2str(i), ' of ' num2str(size(df,1))])

    end
    
    disp('saving.....')
    save('./dataExport/stationSurveyDataTable.mat', 'df')
    
    for w = 1:length(warning_list)
        warning('on', warning_list{w});
    end
    
end

function [data] = calculate_Zu(data, con)
    data.Zu = 0.708 .* ...
        (data.fit.ustar ./ data.DSV) .* ...
        ((data.fit.ustar .* data.grabsumm.d50.*1e-6) ./ con.nu).^0.6;
end

function [modelParams] = XXXXfind_RouseSolns(modelParams, flowDepth, conc)
    modelParams.refNum = (flowDepth - modelParams.b) / modelParams.b; % this is ((H-b)/b), a constant
    modelParams.formula = strcat('sampConcMean ~ cb *( ((', num2str(flowDepth), '-collZ)/collZ) / ', num2str(modelParams.refNum), ' )^Rou');
    
    opts = statset('fitnlm');
    opts.Display = 'off';
    opts.DerivStep = 6.0555e-10;
    options = statset('Display', 'iter', 'TolFun', 1e-10, 'FunValCheck', 'on');
    gof = 0;
    modelParams.guess = [mean(conc.sampConcMean(1:3)) 0.5]; % guesses for cb and Rou
    it = 0;
    while gof < 0.8
        try
            modelParams.model = fitnlm(conc, modelParams.formula, modelParams.guess, 'Options', opts);
            gof = modelParams.model.Rsquared.Ordinary;
        catch
            % try something of a monte carlo approach
            modelParams.guess(1) = mean(conc.sampConcMean(1:3)) + randn(1); % try something of a monte carlo approach
            modelParams.guess(2) = 0.5 + randi(50, 1)/100;
            gof = 0;
        end
    end
    modelParams.Rou = modelParams.model.Coefficients.Estimate(1);
    modelParams.cb = modelParams.model.Coefficients.Estimate(2);
end
