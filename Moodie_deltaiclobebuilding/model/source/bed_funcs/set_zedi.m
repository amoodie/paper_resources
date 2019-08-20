function [zed] = set_zedi(zed_zero_idx, S0, nx, dx)
    %% sets the initial floodplain/delta profile
    % should only be called once at beginning of model initiation
    zed_yi = zed_zero_idx*(S0*dx);
    zed = linspace(zed_yi, zed_yi - S0*(nx*dx), nx+1);
end