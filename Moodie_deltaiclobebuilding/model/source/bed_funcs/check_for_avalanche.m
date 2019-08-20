function [eta] = check_for_avalanche(eta, S, Be, mou_idx, depthick, psi)
	% avalanching condition for stability on the foreset

    [maxS, maxI] = max(S);
    area_move = 1000; % m2
    if maxS > 2e-3
        deta1 = area_move ./ Be(maxI-1);
        deta2 = area_move ./ Be(maxI);
        
        eta(maxI-1) = eta(maxI-1) - deta1;
        eta(maxI) = eta(maxI) + deta2;
    end
    
end