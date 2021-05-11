function [marker_set] = plotWPdata(xdata, ydata, idx, idxOut, fig, varargin)
    %plotWPdata plot Wright Parker data with coloring
    %
    % 
    % just a convenience function for plotting lots of data over and over with slight variations

    colswitch = 'qual';
    [colorSet, col1, col2, col3, col4, col5, col6] = load_colorSet(colswitch);
    mar1 = 'd';
    mar2 = '>';
    mar3 = 'pentagram';
    mar4 = '<';
    mar5 = 'v';
    mar6 = 'hexagram';
    
    
    colors = [col1; col2; col3; col4; col5; col6];

    %% process inputs
    defaultMarker = 's';
    defaultColorVariable = [];
    defaultColorVariableLimits = [];

    p = inputParser;
    addRequired(p, 'xdata');
    addRequired(p, 'ydata');
    addRequired(p, 'idx');
    addRequired(p, 'idxOut');
    addRequired(p, 'fig');
    addParameter(p, 'Marker', defaultMarker);
    addParameter(p, 'ColorByVariable', defaultColorVariable);
    addParameter(p, 'ColorByVariableLimits', defaultColorVariableLimits);
    parse(p, xdata, ydata, idx, idxOut, fig, varargin{:});
    
    %% process varargin here to detemrine how to mark and color
    if any(strcmp('Marker', p.UsingDefaults))
        markers = {mar1; mar2; mar3; mar4; mar5; mar6};
    else
        markers = {p.Results.Marker, p.Results.Marker, p.Results.Marker, ...
                   p.Results.Marker, p.Results.Marker, p.Results.Marker};
    end
    
    if any(strcmp('ColorByVariable', p.UsingDefaults))
        color_flag = 'categorical';
    else
        color_flag = 'variable';
        
        color_var = p.Results.ColorByVariable;
    end

    %% plotting
    idxs = unique(idx);
    marker_set = [];
    % if there are any outliers
    if any(~idxOut)
        for i = 1:length(idxs)
            idxi = ( idx == idxs(i) );
            [~] = plot(xdata(and(idxi, ~idxOut)), ydata(and(idxi, ~idxOut)), ...
                'ks', 'MarkerSize', 6, 'Marker', defMarker, 'LineWidth', 0.5);
        end
    end
    
    for i = 1:length(idxs)
        idxi = ( idx == idxs(i) );
        if strcmp(color_flag, 'categorical')
            marker_set_new = plot(xdata(and(idxi, idxOut)), ydata(and(idxi, idxOut)), ...
                'ks', 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'Marker', markers{i}, 'LineWidth', 0.5);
            marker_set = vertcat(marker_set, marker_set_new);
        elseif strcmp(color_flag, 'variable')
            marker_set_new = plot(NaN, NaN, ...
                'ks', 'MarkerSize', 6, 'Marker', markers{i});
            marker_set = vertcat(marker_set, marker_set_new);
            
            scatter(xdata(and(idxi, idxOut)), ydata(and(idxi, idxOut)), 30, ...
                color_var(and(idxi, idxOut)), 'filled', 'Marker', markers{i}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5)
        end
    end
    
end