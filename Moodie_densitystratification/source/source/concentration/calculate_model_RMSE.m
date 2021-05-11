function [model, RMSE] = calculate_model_RMSE(sample_table, model)
% specific for handling the model structures of the code

    %% total model RMSE
    % get the correct c and z vectors out of the struct
    field_names = fieldnames(model);
    if ismember('CsSum', field_names)
        mcs = model.CsSum;
        mzs = model.ZsSum;
    else
        mcs = model.Cs;
        mzs = model.Zs;
    end
    
    % pull out the data from the table for convenience
    sz = sample_table.sampleZ;
    sc = sample_table.concNW;
    
    % interpolate the model cs at the correct zs for the data
    mc = interp1(mzs, mcs, sz);
    
    total_RMSE = model_RMSE(sc, mc);
    model.RMSE = total_RMSE;
    model.RMSEnorm = total_RMSE ./ mcs(1);
    
    
    %% grain size specific RMSE
    if ismember('Cs', field_names) && ismember('CsSum', field_names)
        % pull out model c_i 
        mcis = model.Cs;
        mzis = model.Zs;
        
        % sample specific c_i
        sc = sample_table.concNWbyClass;
        sc = [sc{:}]';
        
        class_RMSE = NaN(1, size(mcis, 2));
        class_E = NaN(size(sc, 1), size(mcis, 2));
        class_class = [sample_table.gsClass{:}]';
        class_DistNorm = NaN(size(sc, 1), size(mcis, 2));
        
        % loop through each class and calculate RMSE
        for i = 1:size(mcis, 2)
            % interpolate the model cs at the correct zs for the data
            mci = interp1(mzis(:,i), mcis(:,i), sz);
            
            % RMSE calc
            class_RMSE(i) = model_RMSE( sc(:,i), mci );
            class_E(:, i) = (sc(:,i) - mci); % ./ model.params.cb(i);
            class_DistNorm(:, i) = sample_table.sampleZnorm;
        end
        
        model.class_RMSE = class_RMSE;
        model.class_RMSEnorm = class_RMSE ./ mcis(1,:);
        model.class_E = class_E;
        model.class_class = class_class;
        model.class_DistNorm = class_DistNorm;
        
    end

end

function [RMSE] = model_RMSE(sc, mc)
    % calculate RMSE between two
    RMSE = sqrt( mean( (sc - mc) .^2) );
end
