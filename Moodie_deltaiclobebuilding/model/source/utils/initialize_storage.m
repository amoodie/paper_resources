function [s] = initialize_storage(nx, Tdays)

    s.eta = zeros(nx+1, Tdays); % preallocate eta storage matrix
    s.H = zeros(nx+1, Tdays); % preallocate H storage matrix
    s.zed = zeros(nx+1, Tdays); % preallocate zed storage matrix
    s.avul = NaN(1, Tdays);
    s.mou_idx = zeros(1, Tdays);
    s.rad = zeros(1, Tdays);
    s.levee = zeros(nx+1, Tdays);
    s.Lb = NaN(1, Tdays);
    s.dts = zeros(1, Tdays);
    s.Be = zeros(nx+1, Tdays);
    s.Qw = NaN(nx+1, Tdays);

end