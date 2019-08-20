function [Hn] = Hn_DW(Q, B, Cf, S0, g)
	%% generalized Darcy Weisbach, from Ganti et al., 2014 supp info
	Hn = (Q ./ B).^(2/3) .* (Cf./(g.*S0)).^(1/3);

end