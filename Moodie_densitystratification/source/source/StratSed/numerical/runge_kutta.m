function [y]= runge_kutta(f, y0, h)
    
    k1 = f(y0);  % Approx for y gives approx for deriv
    y1 = y0 + k1 * h/2;      % Intermediate value (using k1)

    k2 = f(y1);        % Approx deriv at intermediate value.
    y2 = y1 + k2 * h/2;      % Intermediate value (using k2)

    k3 = f(y2);        % Another approx deriv at intermediate value.
    y3 = y2 + k3 * h;        % Endpoint value (using k3)

    k4 = f(y3);        % Approx deriv at endpoint value.

    y = y0 + (k1 + 2*k2 + 2*k3 + k4) * h/6; % Approx soln at x+h
    
end