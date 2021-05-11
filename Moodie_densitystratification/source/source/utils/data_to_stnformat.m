function [stnform, dist] = data_to_stnformat(data, r)
    nSamp = size(data.nomDepth, 2)-1;
    flowDepth = repmat(data.flowDepth(r), nSamp, 1);
    nomDepth = data.nomDepth(r, 2:end)';
    collDepth = data.collDepth(r, 2:end)';
    collZ = data.collZ(r, 2:end)';
    sampConc = data.sampConc(r, 2:end)';
    dist = data.suspdist{:, r+1};
    if isfield(data, 'suspdistnearonly')
        didxend = (r*3) + 1;
        dist = mean(data.suspdistnearonly{:, didxend-2:didxend}, 2);
    end
    dist(isnan(dist)) = 0;
    stnform = array2table(horzcat(flowDepth, nomDepth, collDepth, collZ, sampConc), ...
        'VariableNames', {'flowDepth', 'nomDepth', 'collDepth', 'collZ', 'sampConc'});
end