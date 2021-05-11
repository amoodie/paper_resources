function [simple] = simpleModel(u_star, C_bed_frac, V, S_0, h, b, Nz)

    %% constants
    g = 9.81; % gravity
    R = 1.65; % specific gravity
    nu = 1.004 * 1e-6; % viscosity
    kappa = 0.41; % von Karman constant

    %% set up z
    z = linspace(b*h, h, Nz)';
    eta = z ./ h;

    %% set up other params
    u_star_1 = u_star / u_star;
    Re_tau = sqrt(g * S_0 * h * h * h) / nu;
    k_s = 0.01; % roughness
    D_frac = length(C_bed_frac);
    C_mu = 0.09;

    %% Without stratification (log and Vanoni)
    Re_star = Re_tau*u_star_1*k_s;
    
    % log
    u_non1 = 1/kappa*   (  log(eta'/k_s) + kappa *  (8.5+(1 / kappa * log(Re_star) - 3) * exp(-0.121 * (log(Re_star))^2.42))  );
    u_non2 = (u_star / kappa) .* log(z ./ k_s);
    u_non2_norm = u_non2 ./ u_star;
    u_non = u_non2_norm';
    
    nu_t_non = kappa*u_star_1*eta'.*(1-eta');
    c_i_non = zeros(Nz,D_frac);

    % Vanoni
    for i = 1:D_frac
        c_i_non(:,i) = C_bed_frac(i)*((1-eta')./eta'*b/(1-b)).^(V(i)/kappa);
    end
    c_non = sum(c_i_non, 2);
    
    %% pack it up
    simple.c_non = c_non;
    simple.c_i_non = c_i_non;
    simple.u_non = u_non;
    simple.nu_t_non = nu_t_non;
    simple.eta = eta;
    simple.z = z;

end