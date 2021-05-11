function [marker_set] = plotYRdata(xdata, ydata, idxYr, idxOut, fig, varargin)
    %plotYRdata plot YR data with outliers and yr coloring
    %
    % 
    % just a convenience function for plotting lots of data over and over with slight variations

    global col1 col2 col3 mar1 mar2 mar3
    colors = [col1; col2; col3];
    figure(fig)
    hold on;
    
    %% process inputs
    defaultMarker = 's';
    defaultColorVariable = [];
    defaultColorCategorical = [];
    defaultColorVariableLimits = [];
    defaultColorEdgeOnly = false;
%     defaultColorIdx = true;
    defaultIgnoreOutliers = false;

    p = inputParser;
    addRequired(p, 'xdata');
    addRequired(p, 'ydata');
    addRequired(p, 'idxYr');
    addRequired(p, 'idxOut');
    addRequired(p, 'fig');
    addParameter(p, 'Marker', defaultMarker);
    addParameter(p, 'ColorByVariable', defaultColorVariable);
    addParameter(p, 'ColorByCategorical', defaultColorCategorical);
    addParameter(p, 'ColorByVariableLimits', defaultColorVariableLimits);
    addParameter(p, 'ColorEdgeOnly', defaultColorEdgeOnly);
    addParameter(p, 'IgnoreOutliers', defaultIgnoreOutliers);
    parse(p, xdata, ydata, idxYr, idxOut, fig, varargin{:});
    
    
    %% process varargin here to detemrine how to mark and color
    if any(strcmp('Marker', p.UsingDefaults))
        markers = {mar1; mar2; mar3};
    else
        markers = {p.Results.Marker, p.Results.Marker, p.Results.Marker};
    end
    
    if and(~isempty(p.Results.ColorByCategorical), ~isempty(p.Results.ColorByVariable))
        error('must not specify both options!')
    end
    
    if and(isempty(p.Results.ColorByCategorical), isempty(p.Results.ColorByVariable))
        color_flag = 'year';
    elseif ~isempty(p.Results.ColorByCategorical)
        color_flag = 'categorical';
        color_var = p.Results.ColorByCategorical;
    elseif ~isempty(p.Results.ColorByVariable)
        color_flag = 'variable';
        color_var = p.Results.ColorByVariable;
    else
        error('why?')
    end
    
    marker_lookup = [2015, 2016, 2018];
    
    %% sanity check of vector lengths
    all_same_len = all([(length(xdata)==length(ydata)), (length(xdata)==length(idxYr)), (length(xdata)==length(idxOut))]);
    if ~all_same_len
        error('vectors not the same length')
    end
    
    %% process out any vector lists to long lists for easier plotting

    %% plotting
    yrIdxs = unique(idxYr);
    marker_set = [];
    % if there are any outliers
    if ~p.Results.IgnoreOutliers
        if any(~idxOut)
            for i = 1:length(yrIdxs)
                idxYri = ( idxYr == yrIdxs(i) );
                markerIdx = find(marker_lookup == yrIdxs(i));
                [~] = plot(xdata(and(idxYri, ~idxOut)), ydata(and(idxYri, ~idxOut)), ...
                     'ks', 'MarkerSize', 6, 'Marker', markers{markerIdx});
            end
        end
    end
    
    % plot non-outliers, color by flag selections
    for i = 1:length(yrIdxs)
        idxYri = ( idxYr == yrIdxs(i) );
        markerIdx = find(marker_lookup == yrIdxs(i));
        if strcmp(color_flag, 'year')
            
            if p.Results.ColorEdgeOnly
                marker_set_new = plot(xdata(and(idxYri, idxOut)), ydata(and(idxYri, idxOut)), ...
                    'ks', 'MarkerFaceColor', 'none', 'MarkerSize', 6, 'Marker', markers{markerIdx}, 'MarkerEdgeColor', colors(markerIdx,:));
            else
                marker_set_new = plot(xdata(and(idxYri, idxOut)), ydata(and(idxYri, idxOut)), ...
                    'ks', 'MarkerFaceColor', colors(markerIdx,:), 'MarkerSize', 6, 'Marker', markers{markerIdx});
            end
            marker_set = vertcat(marker_set, marker_set_new);

        elseif strcmp(color_flag, 'variable')
            marker_set_new = plot(NaN, NaN, ...
                'ks', 'MarkerSize', 6, 'Marker', markers{i});
            marker_set = vertcat(marker_set, marker_set_new);
            scatter(xdata(and(idxYri, idxOut)), ydata(and(idxYri, idxOut)), 28, ...
                color_var(and(idxYri, idxOut)), 'filled', 'Marker', markers{markerIdx}, 'MarkerEdgeColor', 'k')
            
        elseif strcmp(color_flag, 'categorical')
            cats = categories(color_var);
            [~, catsortidx] = sort(cellfun(@str2num, cats));
            cats = cats(catsortidx);
            cmap = parula(length(cats));
            for j = 1:length(cats)
                idxCati = color_var == cats(j);
                if p.Results.ColorEdgeOnly
                    marker_set_new = plot(xdata(and(and(idxYri, idxCati), idxOut)), ydata(and(and(idxYri, idxCati), idxOut)), ...
                        'ks', 'MarkerFaceColor', 'none', 'MarkerSize', 6, 'Marker', markers{markerIdx}, 'MarkerEdgeColor', cmap(j,:));
                else
                    marker_set_new = plot(xdata(and(and(idxYri, idxCati), idxOut)), ydata(and(and(idxYri, idxCati), idxOut)), ...
                        'ks', 'MarkerFaceColor', cmap(j,:), 'MarkerSize', 6, 'Marker', markers{markerIdx});
                end
                if i==1
                    marker_set = vertcat(marker_set, marker_set_new);
                end
                
            end
        end
    end

end