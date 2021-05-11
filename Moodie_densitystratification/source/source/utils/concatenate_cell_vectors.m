function [list] = concatenate_cell_vectors(cell_vectors)
% this function is to grab vectors from structure array and make one long list of the field specified

    max_dims = cellfun(@(x) max(size(x)), cell_vectors);
    if all(max_dims(1) == max_dims(2:end))
        % all dimensions are the same
        shape_vector = size(cell_vectors{1});
        npc = max(shape_vector); % numel per cell

        tp = shape_vector(2) > shape_vector(1); % logical to transpose or not
        list = NaN(npc * size(cell_vectors, 1), 1);
        for i = 1:size(cell_vectors, 1)
            if tp
                list(((i-1)*npc+1):(i*npc)) = cell_vectors{i}';
            else
                list(((i-1)*npc+1):(i*npc)) = cell_vectors{i};
            end
        end
        
    else
        % not all dimensions are the same
        list = NaN(sum(max_dims),1);
        it = 1;
        for i = 1:size(cell_vectors, 1)
            shape_vector = size(cell_vectors{i});
            npc = max(shape_vector); % numel per cell

            list(it:(it+npc-1)) = cell_vectors{i};
            it = it + npc;
        end
        
    end
        
end