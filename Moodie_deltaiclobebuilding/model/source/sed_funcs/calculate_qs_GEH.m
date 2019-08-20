function [qs] = calculate_qs_GEH(U, Cf0, D50, p1, p2, con)
    % backwater flow taua constant Cf
    taua = (Cf0 .* (U.*U)) ./ (con.Rr * con.g * D50);
    qs_n = (p1 .* (taua .^p2)) ./ Cf0; % qs_n == qs*, i.e. non-dimensionalized
    qs = qs_n .* sqrt(con.Rr * con.g * (D50*D50*D50)); % dimensionalize
end