function [modelEvalZs, modelEvalVs] = evaluate_LotwModel(flowDepth, z0, ustar, alpha)
    nEvalPts = 50;
    modelEvalZs = linspace(z0, flowDepth, nEvalPts); % model evaluation points
    modelEvalVs = (ustar / (alpha * 0.41)) .* log(modelEvalZs ./ z0);
end