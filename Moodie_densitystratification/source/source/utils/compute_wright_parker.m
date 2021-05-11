function [alpha] = compute_wright_parker(cb5, S0)

    cb5S = cb5 ./ S0;
    
    alpha = nan(size(cb5S, 1), 1);
    alpha(cb5S <= 10) = 1 - (0.06 .* (cb5S(cb5S <= 10)) .^ 0.77);
    alpha(cb5S > 10) = 0.67 - (0.0025 .* (cb5S(cb5S > 10)));
    
end