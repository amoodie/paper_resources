function [comparison] = plotCompareYears(data15, data16)
    global col1
    global col2
    grabcum15 = cumsum(data15.grabdist{:, 2:end}); % grab cum
    grabcum16 = cumsum(data16.grabdist{2:end, 2:end});
    for i = 1:size(grabcum15, 2)
        wlcut15(i) = find(grabcum15(:, i) > 5, 1, 'first') - 1;
        wlsize15(i) = data15.grabdist{wlcut15(i), 1};
        wldist15 = (data15.suspdist{1:wlcut15(i), i+1});
        wldist15norm = wldist15 ./ sum(wldist15);
        wld5015(i) = interp1(cumsum(wldist15norm(find(wldist15norm > 0, 1):end)), data15.grabdist{find(wldist15norm > 0, 1):wlcut15(i), 1}, 0.5);
        wld9015(i) = interp1(cumsum(wldist15norm(find(wldist15norm > 0, 1):end)), data15.grabdist{find(wldist15norm > 0, 1):wlcut15(i), 1}, 0.9);
    end
    for i = 1:size(grabcum16, 2)
        wlcut16(i) = find(grabcum16(:, i) > 5, 1, 'first') - 1;
        wlsize16(i) = data16.grabdist{wlcut16(i), 1};
        wldist16 = (data16.suspdist{1:wlcut16(i), i+1});
        wldist16norm = wldist16 ./ sum(wldist16, 'omitnan');
        wld5016(i) = interp1(cumsum(wldist16norm(find(wldist16norm > 0, 1):end)), data16.grabdist{find(wldist16norm > 0, 1):wlcut16(i), 1}, 0.5);
        wld9016(i) = interp1(cumsum(wldist16norm(find(wldist16norm > 0, 1):end)), data16.grabdist{find(wldist16norm > 0, 1):wlcut16(i), 1}, 0.9);
    end
    fig = figure();
    plot(repmat(1, size(wld5015',1), 1)+(randn(size(wld5015'))/50), [wld5015'], 'o', 'Color', col1, 'MarkerSize', 4)
    hold on
    plot(repmat(2, size(wld5015',1), 1)+(randn(size(wld5015'))/50), [wld5016'; NaN(19,1)], 'o', 'Color', col2, 'MarkerSize', 4)
    boxplot([wld5015', [wld5016'; NaN(19,1)]])
        ylabel('D_{50} grain size (\mu m)')
        set(gca, 'XTick', 1:2, 'XTickLabel', {'2015', '2016'})
        boxes = findobj(gca, 'tag', 'Box');
        medians = findobj(gca, 'tag', 'Median');
        outliers = findobj(gca, 'tag', 'Outliers');
        set(boxes, 'LineWidth', 1, 'Color', [0 0 0])
        set(medians, 'LineWidth', 1, 'Color', [0 0 0])
        set(outliers, 'LineWidth', 1, 'Color', [0 0 0])
        box on
        set(gca, 'FontSize', 10, 'LineWidth', 1.5)
        pause(0.1)
    set(fig, 'Pos', [100 100 300 250]);
    set(fig, 'PaperPositionMode', 'auto')
    print('-depsc','-r300','-painters', './figs/distcomp_d50susp.eps');
    
    nidx = 200;
    bedidx15 = randi(size(data15.grabdist, 2)-1, nidx, 1)+1;
    bedidx16 = randi(size(data16.grabdist, 2)-1, nidx, 1)+1;
    nearidx15 = randi(size(data15.neardist, 2)-1, nidx, 1)+1;
    nearidx16 = randi(size(data16.suspdistnearonly, 2)-1, nidx, 1)+1;
    bedcmp = NaN(nidx, 2);
    nearcmp = NaN(nidx, 1);
    for i = 1:nidx
        [bedcmp(i, 1), bedcmp(i, 2)] = kstest2(cumsum(data15.grabdist{:, bedidx15(i)}), cumsum(data16.grabdist{2:end, bedidx16(i)}));
        nearcmp(i) = kstest2(cumsum(data15.neardist{:, bedidx15(i)}), cumsum(data16.suspdistnearonly{2:end, nearidx16(i)}));
    end
    comparison = NaN;
    
    t = figure();
    for i = 1:size(data15.neardist, 2)-1
        semilogx(data15.neardist{:, 1}, cumsum(data15.neardist{:, i+1}), 'Color', col1)
        hold on
    end
    for i = 1:size(data16.suspdistnearonly, 2)-1
        semilogx(data16.suspdistnearonly{2:end, 1}, cumsum(data16.suspdistnearonly{2:end, i+1}), 'Color', col2)
    end
    
    beddistcomp = figure();
    for i = 1:size(data15.grabdist, 2)-1
        l1 = semilogx(data15.grabdist{:, 1}, cumsum(data15.grabdist{:, i+1}), 'Color', col1, 'LineWidth', 1.2);
        hold on
    end
    for i = 1:size(data16.grabdist, 2)-1
        l2 = semilogx(data16.grabdist{2:end, 1}, cumsum(data16.grabdist{2:end, i+1}), 'Color', col2, 'LineWidth', 1.2);
    end
    legend([l1 l2], {'2015', '2016'})
    xlim([1e-1 1e3])
    ylim([0 100])
    xlabel('D (\mu m)')
    ylabel('cumulative %')
    box on
    set(gca, 'FontSize', 10, 'LineWidth', 1.5)
    set(beddistcomp, 'Pos', [100 100 600 250]);
    set(beddistcomp, 'PaperPositionMode', 'auto')
    print('-depsc','-r300','-painters', './figs/distcomp.eps');
    
    filename = strcat('./2016_summer_data/', 'yrs_cmp.csv');
    fid = fopen(filename);
    data = csvread(filename, 2, 2);
    data(:, 5) = NaN(size(data(:,5)));
    hdr_raw = textscan(fid, '%s', (size(data, 2)+2), 'delimiter', ',');
    hdr = hdr_raw{:}';
    class_rng(1) = find(strcmp(hdr, 'Result Between User Sizes (Sizes in um)'));
    class_rng(2) = find(strcmp(hdr, 'Result Below User Sizes (Sizes in um)'));
    class_raw = textscan(fid, '%s', (size(data, 2)+2), 'delimiter', ',');
    class = cellfun(@str2double, cell(class_raw{:}(class_rng(1)-1:class_rng(2)-2)')); % EXPLICIT TO NUMBER OF CLASSES!!!
    frewind(fid);
    list_raw = textscan(fid, '%s %s %*[^\n]', 'delimiter', ',');
    list = list_raw{:,2}(3:end);
    distrng = [5 5+length(class)-1];
    
    homocomp = figure();
    semilogx(repmat(class, size(data, 1), 1)', data(:, distrng(1):distrng(2))', 'Color', [0 0 0])
    xlim([1e-1 1e3])
    ylim([0 15])
    xlabel('D (\mu m)')
    ylabel('volume %')
    box on
    set(gca, 'FontSize', 10, 'LineWidth', 1.5)
    set(homocomp, 'Pos', [100 100 350 250]);
    set(homocomp, 'PaperPositionMode', 'auto')
    print('-depsc','-r300','-painters', './figs/homogenizedcomp.eps');
    
end