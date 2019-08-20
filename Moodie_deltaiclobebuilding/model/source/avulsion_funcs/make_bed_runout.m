function [runout, runout_idx] = make_bed_runout(eta, etai, datum, toe_idx, rad_idx)
    % make a bed runout for an avulsion, starts the runout at the toe_idx
    % at the end of the taper input, calculates the runout inside of
    % rad_idx as the higher of datum or etai and etai for outside of
    % rad_idx

    runout_idx = false(size(eta));
    runout_idx(toe_idx:end) = true;
    
    if rad_idx <= toe_idx
        runout = etai(runout_idx);
        
    else
        temp = eta; % just start here
        temp(rad_idx:end) = etai(rad_idx:end); % beyond radius has to be original
        temp(toe_idx:rad_idx-1) = max([etai(toe_idx:rad_idx-1); datum(toe_idx:rad_idx-1)], [], 1); % between toe and radius is bankfull below
        
        runout = temp(runout_idx);
        
    end
end