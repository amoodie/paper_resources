function [velocity] = MakeVelocity18(stnInfo, stationStruct)
    % load in and plot the velocity data
    
    % load velocimeter
    stnInfoSplit = strsplit(stnInfo, '/');
    stnInfoMod = [stnInfoSplit{1} stnInfoSplit{2} stnInfoSplit{3}];
    try
        
        if strcmp(stnInfoMod, 'KC1_070818')
            [meter] = load_meter_swoffer(stnInfoMod); % this station is swoffer style
        else
            [meter] = load_meter_chinese(stnInfoMod);
        end
        
    catch
        [meter] = []; % make empty holder
        warning('load_meter_xxxxx failed in MakeVelocity18 for %s', stnInfo)
    end
    
    % load adcp
    stnInfoSplit = strsplit(stnInfo, '_');
    stnInfoSplitDate = strsplit(stnInfoSplit{2}, '/');
    stnInfoMod = [stnInfoSplit{1} '-' ...
        stnInfoSplitDate{1} stnInfoSplitDate{2} stnInfoSplitDate{3}];
    
    
    % handle various adcp depths...currently breaks due to depth being too
    % small and producing negative measZ values...
    [adcp] = load_adcp(stnInfoMod, 0.42);
    [adcp] = preprocess_adcp(adcp);
    [adcp] = process_adcp(adcp);
    
    % adjust all the meter data to match up with the adcp data:
    if ~isempty(adcp.depth)
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
    
    if isempty(velocity.adcp.depth)
        velocity.adcp.measZ = [];
        velocity.adcp.mean = [];
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
    
    % depth-averaged velocity from measurements
    [sZ, si] = sort(velocity.data.measZ);
    sZm = sZ(1:end-1) + (  (sZ(2:end) - sZ(1:end-1)) / 2  ); % midpoints between data
    sZe = [0; sZm; velocity.data.flowDepth];
    weight = ( sZe(2:end) - sZe(1:end-1) ) ./ velocity.data.flowDepth;
    velocity.data.mean = sum( velocity.data.velMean(si) .* weight ) / sum(weight);

end