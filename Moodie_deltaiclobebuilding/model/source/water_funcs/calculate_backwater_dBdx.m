function [H] = calculate_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx)
    % backwater formulated for changing width
    g = 9.81; % gravitational acceleration constant
    H = NaN(1,nx+1); % preallocate depth 
    H(nx+1) = abs(H0 - eta(nx+1)); % water depth at downstream boundary
    [dBdx] = (B(2:nx+1) - B(1:nx)) ./ dx; % change in width downstream
    for i = nx:-1:1
        Bi = B(i);
        dBdxi = dBdx(i);
        Hip1 = H(i+1);
        % predictor step: computation of a first estimation of the water depth Hp
        Frsqp = ( Qw*Qw / (g * (Bi*Bi) * (Hip1*Hip1*Hip1)) );
        [dHdxp] = ((S(i+1)-(Cf*Frsqp))/(1-Frsqp)) + (Frsqp/(1-Frsqp)) * (H(i+1)/Bi) * (dBdxi);
        Hp = H(i+1) - dHdxp * dx; % solve for H prediction
        % corrector step: computation of H
        Frsqc = ( Qw*Qw / (g * (Bi*Bi) * (Hp*Hp*Hp)) );
        [dHdxc] = ((S(i)-(Cf*Frsqc))/(1-Frsqc)) + (Frsqc/(1-Frsqc)) * (Hp/Bi) * (dBdxi);
        % convolution of prediction and correction, trapezoidal rule
        H(i) = H(i+1) - ( (0.5) * (dHdxp + dHdxc) * dx );
    end
    if Qw < 75
        [~, idx] = find(eta <= 0, 1);
        H(eta >= -0.2) = H(idx + 1);
    end
end