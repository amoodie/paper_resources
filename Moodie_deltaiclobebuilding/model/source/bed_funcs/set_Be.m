function [Be] = set_Be(Bc0, Bf0, Bo0, mou_idx, rad_idx, nx)
	%% set the depositional width over model domain
    
    % preallocate (not needed)
    Be = zeros(1, nx+1);
    
    % up to the delta radius is channel + floodplain
    Be(1:rad_idx) = Bc0 + Bf0;
    
    % lobe width is set to floodplain width behind mouth up to taper
    taper_offset = 0;
    Be(rad_idx:mou_idx-taper_offset) = Bc0 + Bf0; 
    
    % taper from the 4 behind mouth out to the lobe width at mouth
    Be(mou_idx-taper_offset+1:mou_idx) = linspace(Bc0 + Bf0, Bc0 + Bo0, taper_offset);
    
    % beyond the mouth is channel + lobe
    Be(mou_idx+1:end) = Bc0 + Bo0; % constant width beyond mouth
    
end