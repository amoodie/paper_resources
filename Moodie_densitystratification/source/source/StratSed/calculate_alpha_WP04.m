function [WPalpha] = calculate_alpha_WP04(cb, S)

%     if length(cb) ~= length(S)
%         error('mismatched lenghts')
%     end

    cbS = cb ./ S;
    if cbS <= 10
        WPalpha = (1 - (0.06 .* (cbS) .^ 0.77));
    else
        WPalpha = (0.67 - (0.0025 .* (cbS)));
    end
   
    if any(WPalpha < 0)
        error('alpha < 0, probably input cb as g/L')
    end
    
end