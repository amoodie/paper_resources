function plotRouseModelsIndiv(df)
% plot each rouse model into an individual figure
    for p = 1:size(df, 1)
        figure();
        semilogx(df(p).waterSamplesTable.concNW, df(p).waterSamplesTable.sampleZ, 'ro', 'MarkerSize', 5)
        hold on;
        semilogx(df(p).rouseModel.Cs, df(p).rouseModel.Zs, 'k-s', 'MarkerSize', 3)
        semilogx(df(p).rousePred.CsSum, df(p).rousePred.Zs, 'k--^', 'MarkerSize', 3)
        xlabel('conc (g/L)')
        ylabel('dist. above bed (m)')
        legend('samples', 'best-fit', 'pred.')
        box on; set(gca, 'LineWidth', 1.5, 'FontSize', 10)
    end
end