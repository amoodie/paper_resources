function [discharge] = load_discharge(datetime, table)
    
    %% interpolate
    q = interp1(datenum(table.CollectionDatetime), table.discharge, datenum(datetime));
    
    if isnan(q)
        % interpolate across all the nans in 
        table.dischargeCont = table.discharge;
        nanx = isnan(table.dischargeCont);
        t = 1:numel(table.dischargeCont);
        table.dischargeCont(nanx) = interp1(t(~nanx), table.dischargeCont(~nanx), t(nanx));
        
        % reinterpolate now
        q = interp1(datenum(table.CollectionDatetime), table.dischargeCont, datenum(datetime));
    end
    
    discharge = q;
    
        
end