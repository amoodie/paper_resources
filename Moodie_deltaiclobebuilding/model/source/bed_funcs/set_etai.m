function [eta] = set_etai(zed, Hbf, S0, dx)
	%% sets the initial bed profile
	% should only be called once at beginning of model initiation
    eta = zed - Hbf;
    basin = eta < -18;
    Nbasin = length(find(basin));
    brk = eta(find(~basin, 1, 'last'));
    eta(basin) = linspace(brk, brk - (S0*1e-1)*(dx*Nbasin), Nbasin);
end