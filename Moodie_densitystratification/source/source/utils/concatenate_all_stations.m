function [catted] = concatenate_all_stations(df, tocat)


    splt = strsplit(tocat, '.');
    layer = df;
    for i = 1:length(splt)
        if isstruct( eval(['layer(1).', splt{i}]) )
        layer = vertcat(arrayfun(@(x) eval(['x.', splt{i}]), layer));
%         layer = vertcat( eval(['layer.', splt{i}]) )
        else
            catted = NaN(size(layer, 1), size(eval(['layer(1).', splt{i}]), 2));
            for j = 1:size(catted, 1)
                catted(j, :) = eval(['layer(',num2str(j),').', splt{i}]);
            end
        end
        
    end

end