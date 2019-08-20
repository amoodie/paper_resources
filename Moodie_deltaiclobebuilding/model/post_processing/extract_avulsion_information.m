function [d] = extract_avulsion_information(s, d)
    % extract information from the timeseries surrounding the avulsion
    % locations
    
    %% unpack vars if needed
    dx = d.dx;

    d.avul_idx = find(~isnan(s.avul)); % time index of avulsions
    
    d.avul_len = (s.mou_idx(d.avul_idx-1) - s.avul(d.avul_idx))*dx; % length of avulsions
    d.avul_len_mean = mean(d.avul_len)/1000;
    
    d.lobe_len = (s.mou_idx(d.avul_idx-1)*dx - s.rad(d.avul_idx-1)); % length of lobes at avulsion time
    d.lobe_len_mean = mean(d.lobe_len)/1000;
    
    d.avul_time = [NaN (d.avul_idx(2:end) - d.avul_idx(1:end-1))]; % time between avulsions
    d.avul_time_mean = nanmean(d.avul_time(4:end))/365;
    
end