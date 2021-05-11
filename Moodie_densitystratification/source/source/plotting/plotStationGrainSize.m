function plotStationGrainSize(df)
    [uniqueDepth, ~, uniquePlace] = unique(df.waterSamplesTable.sampleDepth);
    cmapRaw = flipud(parula( size(uniqueDepth, 1) ));
    cmap = cmapRaw(uniquePlace, :);
    figure()
    semilogx(gsDistClass2num(df.gsDistBed), cumsum(df.gsDistBed{:,1}, 'omitnan'), 'k-o', 'LineWidth', 1.5, 'MarkerSize', 2)
    hold on
    for p = 1:size(df.waterSamplesTable, 1)
        semilogx(df.waterSamplesTable.gsClass{p}(2:end), cumsum(table2array(df.waterSamplesTable.gsDistNW{p}(2:end,:)), 'omitnan'), ...
            '-o', 'Color', cmap(p, :), 'LineWidth', 1.5, 'MarkerSize', 2)
%         plot( % plot washload cutoff line
    end
    title(strrep(strtrim(df.StationID), '_', ' '));
    xlim([10 500])
    ylim([0 100])
    xlabel('grain size (\mu m)')
    ylabel('% finer')
    box on;
    set(gca, 'LineWidth', 1.5, 'FontSize', 10)
end