function cleanup_colorbar(colb, labels, varargin)

    %% process inputs
    default_title = [];

    p = inputParser;
    addParameter(p, 'title', default_title);
    parse(p, varargin{:});

    % get limits and determine bin spacings and label points
    lims = colb.Limits;
    nClasses = length(labels);
    binedges = linspace(lims(1), lims(2), nClasses+1);
    dbin = [(binedges(2)-lims(1)), (binedges(2:end)-binedges(1:end-1))];
    ticklocs = binedges(1:end-1) + dbin(1:end-1)/2;
    
    % set new ticks and labels
    colb.Ticks = ticklocs;
    colb.TickLabels = cellstr(num2str(labels));
    
    % set up label if given
    colb.Label.Interpreter = 'latex';
    if ~any(strcmp('title', p.UsingDefaults))
        colb.Label.String = p.Results.title;
    end
    
    % set tick interpreter
    colb.TickLabelInterpreter = 'latex';
    
end