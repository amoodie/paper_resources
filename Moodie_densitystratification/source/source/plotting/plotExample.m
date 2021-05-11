function [FIG_example] = plotExample(data15, data16, rvect)
    global col1
    global col2
    stn15.tab = data_to_stnformat(data15, rvect(1));
    stn16.tab = data_to_stnformat(data16, rvect(2));
    stn15.modelEval = horzcat(data15.modelEvalZs(rvect(1), :)', data15.modelEvalCs(rvect(1), :)');
    stn16.modelEval = horzcat(data16.modelEvalZs(rvect(2), :)', data16.modelEvalCs(rvect(2), :)');
    stn15.flowDepth = data15.flowDepth(rvect(1));
    stn16.flowDepth = data16.flowDepth(rvect(2));
    FIG_example = figure();
    subplot(1, 2, 1)
        hold on
        plot(stn15.modelEval(:, 2), stn15.modelEval(:, 1), 'k--', 'LineWidth', 1.5)
        plot(stn15.tab.sampConc, stn15.tab.collZ, 'ok', 'MarkerFaceColor', col1)
        plot(xlim, [stn15.flowDepth stn15.flowDepth], 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 1.5);
            ylim([0 stn15.flowDepth*1.1])
            title('2015 station')
            xlabel('concentration (g/l)')
            ylabel('collection height (m)')
            box on
            set(gca, 'FontSize', 10, 'LineWidth', 1.5)
    subplot(1, 2, 2)
        hold on
        plot(stn16.modelEval(:, 2), stn16.modelEval(:, 1), 'k--', 'LineWidth', 1.5)
        plot(stn16.tab.sampConc, stn16.tab.collZ, 'ok', 'MarkerFaceColor', col2)
        plot(xlim, [stn16.flowDepth stn16.flowDepth], 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 1.5);
            ylim([0 stn16.flowDepth*1.1])
            title('2016 station')
            xlabel('concentration (g/l)')
            ylabel('collection height (m)')
            box on
            set(gca, 'FontSize', 10, 'LineWidth', 1.5)
    set(FIG_example,'Visible', 'on');
    set(FIG_example, 'Pos', [100 100 650 300]);
    set(FIG_example, 'PaperPositionMode', 'auto')
    print('-dpng', '-r300', './figs/Rouse_example.png');
    print('-depsc', '-painters', '-r300', './figs/Rouse_example.eps');
end