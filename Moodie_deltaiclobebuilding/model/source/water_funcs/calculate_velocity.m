function [U] = calculate_velocity(Qw, H, B)
    U = Qw ./ (H .* B);
end