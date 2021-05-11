function [modelEvalZsN, modelEvalVsN] = normalize_model(modelEvalZs, modelEvalVs)
    % generic normalization by the bottom value, Vs is an arbitrary Value
    
    modelEvalZsN = modelEvalZs ./ modelEvalZs(end);
    modelEvalVsN = modelEvalVs ./ modelEvalVs(1);
    
end