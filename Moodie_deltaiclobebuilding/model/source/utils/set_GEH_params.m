function [p1, p2] = set_GEH_params(def)
    switch def
        case 'EH'
        % use original Engelund and Hansen, 1967 relation for bed matrial load
            p1 = 0.05; % EH original alpha
            p2 = 2.5; % EH original n
        case 'HM'
        % use modfied EH from Ma et al., 201? relation for YR
            p1 = 0.895; % HMqt alpha, Lijin
            p2 = 1.678; % HMqt n, Lijin
        case 'adj'
        % set to whatever, just exploring parameter space here
            p1 = 0.5; % adjust for reasonable
            p2 = 2; % adjust for reasonable
    end
end
