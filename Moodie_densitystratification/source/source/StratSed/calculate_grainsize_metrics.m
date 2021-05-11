function [gsMetrics] = calculate_grainsize_metrics(u_star, D_dim, C_bed_frac)
    
    %% constants
    g = 9.81; % gravity
    R = 1.65; % specific gravity
    nu = 1.004 * 1e-6; % viscosity
    kappa = 0.41; % von Karman constant
    con.g = g; % gravitational constant
    con.rho_f = 1000; % fluid density, kg/m^3
    con.rho_s = 2650; % particle density, kg/m^3
    con.nu = nu; % fluid kinematic viscosity, m^2/s
    
    D_frac = length(C_bed_frac);
    C_frac = C_bed_frac;

    %% metrics
    Re_p = (R * g * D_dim).^0.5 .* D_dim / nu;    
    D_g_bed = 2^(sum(log2(D_dim) .* C_bed_frac));
    sigma_g_bed = 2^( (sum(  ( (log2(D_dim) - log2(D_g_bed)).^2 ) .* C_bed_frac  )) ^0.5 );
    lambda_str = 1 - 0.288 * log2(sigma_g_bed); %WP04 lambda
    
    v_s = ( get_DSV(D_dim, 0.7, 3.5, con) );
    V =  v_s(1:D_frac) ./ u_star; % dimensionless settling velocity
    v_hat = V;
    D_finer(1,1) = 10^(-5);
    per_finer(1,1) = 0;
    if length(D_dim) > 1
        for i = 1:D_frac
            D(i) = fzero(@(D2) sol_D_star(g, R, nu, v_s(i), D2), [10^(-6), 1000*10^(-6)]);
            D_finer(i+1) = (D(i))^2/D_finer(i);
            per_finer(i+1) = per_finer(i)+C_frac(i)*100;
        end
    else
        D = 0;
          
    end
    
    D_g = 2^(sum(log2(D).*(C_frac)));
    sigma_g = 2^( (  sum( ( ( log2(D) - log2(D_g) ).^2) .* (C_frac) )  )^0.5 );
    
    A_ent = 1.30E-07;
    m_hiding = 0.2;
    Zu = u_star./(V).*Re_p.^0.6.*(D/D_g_bed).^m_hiding;
    Ent = A_ent*(lambda_str*Zu).^5./(1+A_ent/0.3*(lambda_str*Zu).^5);

    %% pack into struct
    gsMetrics.Re_p = Re_p;
    gsMetrics.D_g_bed = D_g_bed;
    gsMetrics.sigma_g_bed = sigma_g_bed;
    gsMetrics.lambda_str = lambda_str;
    gsMetrics.v_s = v_s;
    gsMetrics.V = V;
    gsMetrics.D_finer = D_finer;
    gsMetrics.D = D;
    gsMetrics.D_g = D_g;
    gsMetrics.sigma_g = sigma_g;
    gsMetrics.Zu = Zu;
    gsMetrics.Ent = Ent;
    
    
end