function chi_mapping()
    % Andrew J Moodie
    % Oct 2016
    % chi mapping the Appalachians and Gibraltar Arc
    % objects are placed into Matlab structure for batch processing locations
    
    close all;
    APPDEM = GRIDobj('appalachianDEM.tif');
    
    % process the DEM
    [APP] = process_DEM(APPDEM); 
    APP.location = 'Appalachian';
    APP.filename = 'App_'; % filename prefix (for output)
    APP.zone = 17; % UTM zone
    
    % extract streams
    APP.minarea = 200; % minimum thresholding area for channel network
    APP.minlen = 10000; % min channel length
    APP.nKeep = 20; % number of networks to retain (keeps just the largest N)
    [APP.S] = extract_streams(APP);
    
    % calculate chi
    APP.a0 = 1e6; % chi ref area
    APP.mn = 0.45; % chi m/n (theta) ratio, concavity, constant for all networks
    [APP.CHI] = find_chi(APP);
    
    % aggregate data along stream segments and into a mapstruct
    APP.seglen = 20000; % segment length
    [APP.MS] = make_mapstruct(APP);
    
    % load aux data
    fid = fopen('./actual_divide.xy'); frewind(fid); raw = textscan(fid, '%f %f');
    APP.div_act = horzcat(raw{1}, raw{2});
    fid = fopen('./synthetic_divide.xy'); frewind(fid); raw = textscan(fid, '%f %f');
    APP.div_syn = horzcat(raw{1}, raw{2});
    
    % make the map
    APP.binset = 'raw'; % how to bin the data visually
    APP.limsWGS = [-87, 41.5; -74, 33]; % upper left, lower right
    APP.limsUTM = ll2utm(APP.limsWGS(:, 2), APP.limsWGS(:, 1), APP.zone);
    APP.xlim = sort(APP.limsUTM(:, 1))';
    APP.ylim = sort(APP.limsUTM(:, 2))';
    [APP] = make_map(APP);
end


function [OUT] = process_DEM(DEM)
    DEM = fillsinks(DEM);
    DEM.Z(DEM.Z<0) = NaN;      % replace bathymetry with NaN
    
    FD = FLOWobj(DEM, 'preprocess', 'carve');   % flow direction GRIDobj
    A = flowacc(FD);                            % flow accumulation GRIDobj
    G = gradient8(DEM);                         % gradient GRIDobj
    
    DEM = imposemin(FD, DEM, 1e-5); % impose a minimum gradient in the stream network (carving)
    
    % attach to location object for output
    OUT.DEM = DEM;
    OUT.FD = FD;
    OUT.A = A;
    OUT.G = G;
end

function [S] = extract_streams(DAT)
    % modify the stream network to work well with chi analysis
    S = STREAMobj(DAT.FD, 'minarea', DAT.minarea); % find streams
    S = removeshortstreams(S, DAT.minlen); % remove small streams
    S = klargestconncomps(S, DAT.nKeep); % keep largest N
end

function [CHI] = find_chi(DAT)
    CHI = GRIDobj(DAT.DEM, 'double'); % preallocate
    wbar1 = waitbar(0, sprintf('calculating chi of stream %i of %i', 1, DAT.nKeep));
    for s = 1:DAT.nKeep                     % loop through each stream network individually
        T = extractconncomps(DAT.S, s);     % isolate the network T from total network S
        
        chi = chitransform(T, DAT.A, 'mn', DAT.mn); % calculate the chi values following the integral method
        % note that the above calculates chi for the entire DEM
        CHI.Z(T.IXgrid) = chi; % replace stream network with applicable chi values
        
        % update the waitbar
        wbarupdate = s/DAT.nKeep;
        if s < DAT.nKeep
            waitbar(wbarupdate, wbar1, sprintf('calculating chi of stream %i of %i', s+1, DAT.nKeep));
        end
    end
    close(wbar1);
end

function [MS] = make_mapstruct(DAT)
    % create mapstruct of data, aggregating data along a segment
    MS = STREAMobj2mapstruct(DAT.S, 'seglength', DAT.seglen, 'attributes',...
            {'chi', DAT.CHI @mean ...
            'uparea' (DAT.A .* (DAT.A.cellsize^2)) @mean ...
            'gradient' DAT.G @mean});
end

function [DAT] = make_map(DAT)
    % do some creative Matlab plotting because I don't have the mapping toolbox...
    DAT.chiFIG = figure(); % chi mapping figure
    disp('plotting........this may take a moment.......')
    imageschs(DAT.DEM, DAT.DEM, 'colormap', gray, 'colorbar', false);
    hold on
    
    % bin each segment into a percentile and color by percentile
    switch DAT.binset
        case 'quantile'
            DAT.bins = quantile(cell2mat({DAT.MS.chi}'), 0.1:.1:1); % find quantiles
        case 'raw'
            DAT.bins = linspace(0, max(cell2mat({DAT.MS.chi}')), 11); % 10 bin evenly by value
    end
    Cmap = jet(10); % make Cmap for chi data
    for i = 1:length(DAT.MS) % loop through each line segment
        col = Cmap(find(DAT.MS(i).chi <= DAT.bins, 1, 'first')-1, :); % find the appropriate color by percentile
        hold on
        plot(DAT.MS(i).X, DAT.MS(i).Y, 'LineWidth', 1, 'Color', col) % plot the line segment
    end
    
    [tx, ty] = ll2utm(DAT.div_act(:,2), DAT.div_act(:,1), DAT.zone);
    plot(tx, ty, 'LineWidth', 2, 'Color', [1 1 1], 'LineStyle', '-')
    
    [tx, ty] = ll2utm(DAT.div_syn(:,2), DAT.div_syn(:,1), DAT.zone);
    plot(tx, ty, 'LineWidth', 2, 'Color', [1 1 1], 'LineStyle', ':')

    xlim(DAT.xlim)
    ylim(DAT.ylim)
    
    % pretty it up and save
    title([DAT.location, ' chi map'])
    xlabel('easting')
    ylabel('northing')
    colormap(jet(10));
    caxis([0 max(cell2mat({DAT.MS.chi}'))])
    cbar = colorbar;
    set(cbar, 'YTick', [DAT.bins], 'YTickLabel', num2cell(round(DAT.bins./1e4)))
    title(cbar, '\chi value \times 10^4')
    box on
    set(gca, 'FontSize', 10, 'LineWidth', 1.5)
    set(cbar, 'FontSize', 10, 'LineWidth', 1.5)
    set(DAT.chiFIG, 'Pos', [100 100 800 600]);
    set(DAT.chiFIG, 'PaperPositionMode', 'auto')
    print('-dpng', '-r300', strcat('./',DAT.filename,'chimap.png'));
    print('-depsc','-r300','-painters', strcat('./',DAT.filename,'chimap.eps'));
end
