function [] = deviation_calculation()
%     this script uses a "noisy" sine curve and a smoother sine curve to
%     simulate the actual and synthetic divide, respectively.

    X = linspace(-5, 2, 200)'; % longitude
    actDiv = [sin(0.6 .* X) + 0.5*sin(3 .* X) + (randi(100, size(X))/75) + (0.5 .* X), X];
    synDiv = [2.*sin(0.7 .* X+80) + (1 .* X) + 0, X];
    
    [alongDev, alongDist, ~, ~] = get_dists(actDiv, synDiv, 1);
    
    sf = 2; % sig figs
    RMSD = get_RMSD(alongDev, sf);
    
    % map the "divides" and the deviation
    figure()
    subplot(3, 1, [1 2])
        hold on
        plot(actDiv(:, 2), actDiv(:, 1), 'k-', 'LineWidth', 1.5)
        plot(synDiv(:, 2), synDiv(:, 1), 'b-', 'LineWidth', 1.5, 'LineStyle', '--')
        xlim([-6 3]); ylim([-10 10]);
        xlabel('lon'); ylabel('lat');
        box on; set(gca, 'FontSize', 10, 'LineWidth', 1.5);
    subplot(3, 1, 3)
        hold on
        plot(alongDist, alongDev)
        plot([0 max(alongDist)], [RMSD(1) RMSD(1)], 'k--')
        xlabel('along divide distance (km)'); ylabel('dev. distance');
        box on; set(gca, 'FontSize', 10, 'LineWidth', 1.5);
        
end

function [gminds, along, lon, lat] = get_dists(act, syn, dir)
    % dir is boundary direction, 1=N-S, 2=E-W
    [lims, chkdir] = get_lims(syn, dir);
    gminds = NaN(1, size(act,1));
    along = NaN(1, size(act,1));
    lon = NaN(1, size(act,1));
    lat = NaN(1, size(act,1));
    along(1) = 0;
    gminds(1) = min(search(act(1,:), syn));
    for i = 2:size(act,1) % for each point in actual
        delX = get_geodist(act(i,:), act(i-1,:)); % get distance from previous point to this point
        along(i) = along(i-1) + delX; % distance along is this plus the previous
        lon(i) = act(i, 1); % save lon of each pt
        lat(i) = act(i, 2); % save lat of each pt
        if act(i,chkdir) < max(lims) && act(i,chkdir) > min(lims); % if the actual point is within the limits of the synthetic
            [gdists] = search(act(i,:), syn); % search along the synthetic for shortest
            gminds(i) = min(gdists); % the minimum of the vector is the shortest
        else
            gminds(i) = NaN; % if false reset to zero
        end 
    end
end

function [lims, chkdir] = get_lims(syn, dir)
    if dir == 1
        lims = [min(syn(:,2)) max(syn(:,2))];
        chkdir = 2;
    elseif dir == 2
        lims = [min(syn(:,1)) max(syn(:,1))];
        chkdir = 1;
    elseif dir == 3
        lims = [NaN NaN];
        chkdir = 1; % NaN?
    else
        error('directional boundary not set correctly')
    end
end

function [gdists] = search(actp, syn)
    gdists = NaN(1, size(syn,1));
    for j = 1:size(syn,1) % for each point in the synthetic divide
        gdists(j) = get_geodist(actp, syn(j,:));
    end
end

function [dist] = get_geodist(P1, P2)
    % takes two x-y pairs: P1 and P2 (x = column 1)
    % in the case of distance from divides, P1 is actual
    % this returns distance as km
    phi1 = deg2rad(P1(2));
    phi2 = deg2rad(P2(2));
    dphi = deg2rad(abs(P2(2) - P1(2)));
    dlam = deg2rad(abs(P2(1) - P1(1)));
    a = sin(dphi/2) * sin(dphi/2) + ... % haversine formula
        cos(phi1) * cos(phi2) * ...
        sin(dlam/2) * sin(dlam/2);
    dsig = 2*asin(sqrt(a)); % haversine formula
    r = 6371; % km radius of earth
    dist = r*dsig;
end

function [RMSD] = get_RMSD(dists, sf)
    RMSDt = sqrt(  sum( dists .^ 2, 'omitnan') / sum(~isnan(dists))  );
    uncertdists = dists;
    perr = 1.41; % position error input for the calculations
    uncertdists(uncertdists < perr) = perr;
    RMSDuncert = sqrt(  sum( (perr ./ uncertdists) .^ 2, 'omitnan') / sum(~isnan(dists))  );
    RMSD = [round(RMSDt, sf) round(RMSDuncert, sf)]; % RMSD +- uncertainty
end