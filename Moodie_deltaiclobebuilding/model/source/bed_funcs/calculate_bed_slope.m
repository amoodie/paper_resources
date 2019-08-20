function [S] = calculate_bed_slope(eta, nx, dx)
    % return slope of input bed (eta)
    S = zeros(1,nx+1);
    S(1) = (eta(1) - eta(2)) / dx;
    S(2:nx) = (eta(1:nx-1) - eta(3:nx+1)) ./ (2*dx);
    S(nx + 1) = (eta(nx) - eta(nx + 1)) / dx;
end