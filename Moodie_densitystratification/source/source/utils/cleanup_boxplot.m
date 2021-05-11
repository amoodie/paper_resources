function cleanup_boxplot(varargin)
    
    %% process inputs
    default_color = [];

    p = inputParser;
    addParameter(p, 'color', default_color);
    parse(p, varargin{:});
    
    % drawnow to update
    drawnow; pause(0.5)
    
    %% get count of objects and sanity check things
    boxes = findobj(gca, 'tag', 'Box');
    nbox = length(boxes);
    if size(p.Results.color, 1) > 1
        % if more than one row, must match on length (or be longer?)
        if size(p.Results.color, 1) ~= nbox
            error('color must have same length as nbox if M>1')
        end
    end
    
    if ~isempty(p.Results.color) && ~isscalar(p.Results.color)
        if ~(size(p.Results.color, 2) == 3)
            error('must be an Mx3 matrix or scalar')
        end
    end
   
    
    %% process varargin here to detemrine how to mark and color
    if any(strcmp('color', p.UsingDefaults))
        colorlist = [0 0 0] .* ones(nbox, 3); % set it to black
    else
        if isscalar(p.Results.color)
            % assume it is a grayscale
            colorlist = repmat(p.Results.color, nbox, 3) .* ones(nbox, 3);
        elseif size(p.Results.color, 1) == 1
            % assume it is an RGB triplet, but needs to be repeated to fill out list
            colorlist = repmat(p.Results.color, nbox, 3) .* ones(nbox, 3);
        else
            % assume it is a valid colormap or list
            colorlist = flipud(p.Results.color); % need to flip because it gets the boxes in reverse order!!
        end
    end
    
    
    % apply
    boxes = findobj(gca, 'tag', 'Box');
    medians = findobj(gca, 'tag', 'Median');
    outliers = findobj(gca, 'tag', 'Outliers');
    set(boxes, 'LineWidth', 1, {'color'}, num2cell(colorlist, 2));
    set(medians, 'LineWidth', 1, {'color'}, num2cell(colorlist, 2));
    set(outliers, 'Marker', 'none')
    box on
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'layer', 'top', 'TickLabelInterpreter', 'latex')     
    
end