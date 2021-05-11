function [modelEvalZs, modelEvalCs] = evaluate_RouseModel(flowDepth, cb, Hbb, Rou)
    nEvalPts = 51;
    modelEvalZs = linspace(flowDepth*0.05, flowDepth, nEvalPts)';
    modelEvalCs = cb .* ( ((flowDepth-modelEvalZs)./modelEvalZs) ./ Hbb ) .^ Rou;
end