function [Bc] = set_Bc(mou_idx, Bc0, thet, nx, dx)
	%% set the hydraulic width over the model domain
    Bc = zeros(1, nx+1); % preallocate needed?
    
    % offset the spreading by a few nodes as a "jet" beyond mouth
    offset = 2;
    jet_idx = mou_idx+offset;
    
    Bc(1:jet_idx) = Bc0;
    Bc(jet_idx+1:end) = Bc0 + 2 * (linspace(1, (nx+1-jet_idx)*dx, (nx+1-jet_idx)) * tand(thet));
    
    % even smoothing across mouth interval + offset
    n = 10;
    n2 = n*2;
    Bc(jet_idx-n:jet_idx+n+1) = conv(Bc(jet_idx-n2:jet_idx+n2), ones(n2, 1) / n2, 'valid');
    
    
end