function [d] = extract_additional_information(s, d)
    % this would be the place to do additional calculations.
    
    
    
    
    
    % NONE of the below code has been reviewed since major changes and probably
    % doesn't work. It is just here for the sake of being here. This
    % function does nothing unless you uncomment the below calls.

    % [d] = calc_backwaterzone(s, d); 
    % [d.agg_lij] = calc_agg(s, d.lij_idx);
    % [d.prog] = calc_prog(s);

end


function [data] = calc_agg(s, idx)
    data.detadt = (s.eta(idx, 2:end) - s.eta(idx, 1:end-1)) * 1000; % mm/d
    data.dzeddt = (s.zed(idx, 2:end) - s.zed(idx, 1:end-1)) * 1000;
    data.thick = s.eta(idx, 2:end) - s.eta(idx, 1);
end

function [rates] = calc_prog(s)
    times = find(~isnan(s.avul));
    dtime = [NaN (times(2:end) - times(1:end-1))];
    rads = s.rad(times);
    rates.prog_avul = rads ./ dtime;
end

function [d] = calc_backwaterzone(s, d)
    nAvul = length(d.avulidx);
    if nAvul == 0
       d.setidx = true(s.Tdays, 1);
    elseif nAvul > 0 && nAvul < 3
        d.setidx = false(s.Tdays, 1);
        d.setidx(1:d.avulidx) = true;
    else
       d.setidx = false(s.Tdays, 1);
       d.setidx(d.avulidx(end-2):d.avulidx(end-1)-1) = true;
    end
    d.setlen = length(find(d.setidx));
    d.cumChange = s.eta(:, find(d.setidx, 1, 'last')-1) - s.eta(:, find(d.setidx, 1, 'first')+1);
%     d.Qwsplit = ((max(s.Qw(d.setidx)) - min(s.Qw(d.setidx))) / 2) + min(s.Qw(d.setidx));
    d.Qwsplit = 1000; % could try playing with this for different results?
    d.lowQws = padarray(s.Qw(d.setidx) < d.Qwsplit, [0 find(d.setidx, 1, 'first')], 'pre');
    d.highQws = padarray(s.Qw(d.setidx) >= d.Qwsplit, [0 find(d.setidx, 1, 'first')], 'pre');
    d.lowStart = ~d.highQws(2:end) & d.highQws(1:end-1);
    d.highStart = ~d.lowQws(2:end) & d.lowQws(1:end-1);

    d.highchg = s.eta(:, d.highQws) - s.eta(:, [d.highQws(2:end) false]);
    d.highchg = d.highchg(:, 1:end-1);
    [~, d.highdepfrt] = max(d.highchg);
    
    d.lowchg = s.eta(:, d.lowQws) - s.eta(:, [d.lowQws(2:end) false]);
    d.lowchg = d.lowchg(:, 1:end-1);
    [~, d.lowdepfrt] = max(d.lowchg);
    
    d.endavul = s.avul(find(d.setidx, 1, 'last')+1);
    d.startavul = s.avul(find(d.setidx, 1, 'first'));
    d.endmou = s.mou(find(d.setidx, 1, 'last'));
    d.startmou = s.mou(find(d.setidx, 1, 'first'));
    
    switchlist = [find(d.setidx, 1, 'first'), find(d.lowStart), find(d.highStart), find(d.setidx, 1, 'last')];
    lowN = length(find(d.lowStart));
    highN = length(find(d.highStart));
    [~, sortidx] = sort(switchlist);
    d.change = s.eta(:, switchlist(sortidx(2:end))) - s.eta(:, switchlist(sortidx(1:end-1)));
    if s.Qw(find(d.setidx, 1, 'first')) < d.Qwsplit % if lowflow start
        d.lowChange = d.change(:, 1:2:size(d.change, 2));
        d.highChange = d.change(:, 2:2:size(d.change, 2));
    else % if highflow start
        d.lowChange = d.change(:, 2:2:size(d.change, 2));
        d.highChange = d.change(:, 1:2:size(d.change, 2));
    end

    d.highMean = mean(d.highChange(:, 2:end-1), 2);
    highStd = std(d.highChange(:, 2:end-1), 0, 2);
    highTvals = [d.highMean + highStd, d.highMean - highStd];
    d.highShade = [max(highTvals, [], 2)' fliplr(min(highTvals, [], 2)')];
    d.lowMean = mean(d.lowChange(:, 2:end-1), 2);
    lowStd = std(d.lowChange(:, 2:end-1), 0, 2);
    lowTvals = [d.lowMean + lowStd, d.lowMean - lowStd];
    d.lowShade = [max(lowTvals, [], 2)' fliplr(min(lowTvals, [], 2)')];
end

