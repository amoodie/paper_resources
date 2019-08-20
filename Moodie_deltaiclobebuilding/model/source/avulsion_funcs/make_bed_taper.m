function [taper, taper_idx] = make_bed_taper(eta, datum, loc, win, dx)
    % make a bed taper for an avulsion, centers the taper on loc at hwin
    % distance upstream and downstream. Tapers down to the level given by
    % loc+hwin in vector datum

    hwin = floor((win-1)/2);
    taper = NaN(1, win);
    taper_idx = false(size(eta));
    try
        m = (eta(loc-hwin) - datum(loc+hwin)) / (dx*win);
        for i = 1:win
            taper(i) = eta(loc-hwin) - (m * dx * i);
        end
    catch
        taper = eta(loc-hwin:loc+hwin);
    end
    taper_idx(loc-hwin:loc+hwin) = true;
end