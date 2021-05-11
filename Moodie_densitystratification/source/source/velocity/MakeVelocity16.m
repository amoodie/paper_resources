function [velocity] = MakeVelocity16(stnInfo, stationStruct)
    % load in and plot the velocity data
    
    % load velocimeter
    stnInfoSplit = strsplit(stnInfo, '/');
    stnInfoMod = [stnInfoSplit{1} stnInfoSplit{2} stnInfoSplit{3}];
    try
        [meter] = load_meter_swoffer(stnInfoMod);
    catch
        [meter] = []; % make empty holder
    end
    
    % load adcp
    stnInfoSplit = strsplit(stnInfo, '_');
    stnInfoMod = [stnInfoSplit{1}];
    tddepth = 0.4;
    [adcp] = load_adcp(stnInfoMod, tddepth);
    [adcp] = preprocess_adcp(adcp);
    [adcp] = process_adcp(adcp);
   
    % adjust all the meter data to match up with the adcp data:
    %   the adcp always had larger values, we'll go with this since the meter was 
    %   on a pole that was bent through the flow depth.
    if ~isempty(meter)
        minAdcp = stationStruct.FlowDepthAtCollection - max(adcp.depth);
        maxMeter = max(meter.measZ);
        meterMatch = find(meter.measZ > minAdcp);
        adcpMatch = find((stationStruct.FlowDepthAtCollection - adcp.depth) < maxMeter);
        if and(~isempty(meterMatch), ~isempty(adcpMatch))
            % just add a constant over the whole depth for now, could change to a function later
            meterVals = meter.velmean(meterMatch);
            adcpVals = adcp.mean(adcpMatch);
            adjVal = mean(adcpVals) - mean(meterVals);
            
            % apply the adjustment
            meter.velmean = meter.velmean + adjVal;
            meter.velmedian = meter.velmedian + adjVal;
        end
    end
    
    % process data
    [velocity] = velocity_to_struct(meter, adcp);
    
    % if velocity profile is empty, fill with some placeholders
    if isempty(velocity.meter)
        velocity.meter.measZ = [];
        velocity.meter.velmean = [];
    end

    % if adcp profile is empty, fill with some placeholders
    if isempty(velocity.adcp.depth)
        velocity.adcp.measZ = [];
        velocity.adcp.mean = [];
    % otherwise, compute the values to fill out the adcp dataset
    else
        velocity.adcp.flowDepth = stationStruct.FlowDepthAtCollection;
        velocity.adcp.measZ = velocity.adcp.flowDepth - velocity.adcp.depth;
        velocity.adcpTable.measZ = velocity.adcp.measZ';
    end
    
    % collapse to a single velocity dataset for analysis
    velocity.data.measZ = [velocity.adcp.measZ'; velocity.meter.measZ'];
    velocity.data.velMean = [velocity.adcp.mean'; velocity.meter.velmean'];
    velocity.data.source = [repmat({'adcp'}, length(velocity.adcp.mean),1); repmat({'meter'}, length(velocity.meter.velmean),1)];
    velocity.data.flowDepth = stationStruct.FlowDepthAtCollection;
    
    [sZ, si] = sort(velocity.data.measZ);
    sZ = [sZ; velocity.data.flowDepth];
    weight = ( sZ(2:end) - sZ(1:end-1) ) ./ velocity.data.flowDepth;
    velocity.data.mean = sum( velocity.data.velMean(si) .* weight ) / sum(weight);  % depth averaged
end