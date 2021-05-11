function [data] = velocity_to_struct(meter, adcp)
    data.meter = meter;
    data.adcp = adcp;
    if ~isempty(meter)
        data.meterTable = array2table([data.meter.measZ' data.meter.velmean' ...
            data.meter.velmedian'], ...
            'VariableNames', {'measZ', 'mean', 'median'});
        data.meterTable.err = data.meter.velerr';
    end
    if ~isempty(adcp.depth)
        data.adcpTable = array2table([data.adcp.depth', data.adcp.mean' data.adcp.median' data.adcp.std'], ...
            'VariableNames', {'measDepth', 'mean', 'median', 'std'});
    end
end