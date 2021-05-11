function [strat] = stratModel(u_star, D_dim, cb_i, C_bed_frac, V, S_0, b, h, Nz, solution_plot, final_plot, verbose, simple)

    %% set up constants
    sigma_sc = 1; % schmidt number
    
    A_1 = 0.92; % simple constant
    B_1 = 16.6; % simple constant
    y_1 = 0.22;
    C_1 = 0.08; % simple constant
    A_2 = A_1 * (y_1 - C_1) / y_1 / sigma_sc; % works out to nearly same value as WP04
    B_2 = 10.1; % simple constant
    X_1 = A_2 - ((6*A_1*A_2) / B_1); % combinatory constant
    X_2 = (18*A_1*A_2) + (3*A_2*B_2); % combinatory constant, printed as X_2 = (18*A1*As) + (3*A2*B2)
    X_3 = (9*A_1*A_2); % combinatory constant
    X_4 = A_1 - ((6*A_1^2)/B_1) - (3*A_1*C_1); % combinatory constant
    X_5 = (18*A_1^2) + (9*A_1*A_2);
    
    con.g = 9.81; % gravitational constant
    con.rho_f = 1000; % fluid density, kg/m^3
    con.rho_s = 2650; % particle density, kg/m^3
    con.nu = 1.004 * 1e-6; % fluid kinematic viscosity, m^2/s
    g = con.g; % gravity
    R = 1.65; % specific gravity
    nu = con.nu; % viscosity
    kappa = 0.41; % von Karman constant
    
    C_mu = 0.09; % tke scaling constant
    alpha = 0.7; % turbulence length scale constant 
    MO_a = 0; %  monin-obukhov length?
    
    dt_ini = 0.005; % timestep for adv-diff eqn
    dt = dt_ini;
    u_star_1 = 1;
    
    k_s = 0.01; % roughness
    % b = b; % roughness (passed in as arg)
    
    tol = 1e-6; % solution tolerance
    maxiter = 10000; % max number of iterations in solution
    
    %% set up z
    z = linspace(b*h, h, Nz)';
    dz = z(2)-z(1);
    eta = z ./ h;
    deta = eta(2) - eta(1);
    err_idx = (1:Nz) < (floor(Nz/2));
    
    %% calculate any new parameters needed
    D_frac = length(D_dim);
    C_frac = C_bed_frac;
    
    Re_p = (R * g * D_dim).^0.5 .* D_dim / nu; 
    Ri_tau = R * mean([0, sum(cb_i)]) / 2 / S_0; % function rearrangef from given: S = R*C_mean/Ri_tau
    Re_tau = sqrt(g * S_0 * h * h * h) / nu; % function rearranged from given: h = (((Re_tau*nu)^2)/g/S)^(1/3)

    if (Re_tau * k_s > 200)
        HydrRough = 'r';
    elseif (Re_tau * k_s < 3)
        HydrRough = 's';
    else
        HydrRough = 't';
    end
    
    Re_flag = 0;
    Ri_flag = 1;
    f = 1;
    
    one_mat = ones(1, Nz);
    zero_mat = zeros(1, Nz);
    
    parkat = 500;
    scheme = 6;
    beta_TVD = 0;
    
    
    %% With stratification
    % do some preallocation of vectors and then fill with initials
    u_hat = (((u_star / kappa) .* log(z ./ k_s) ) ./ u_star);
    
    c = repmat(sum(cb_i), Nz, 1);
    C_mean = mean(c);
    c_i = repmat(cb_i, Nz, 1) ; % gs specific profiles, volumetric conc
    if C_mean > 0
        c_hat = c ./ C_mean; % total concentration at eta, nondim to cb
        c_hat_i = c_i ./ C_mean; % gs specific profiles, nondim to cb
    else
        c_hat = c;
        c_hat_i = c_i;
    end
    k_hat = (one_mat'.^2) ./ (C_mu .^ 0.5) .* (1 - b); % (  (1.^2) ./ (C_mu .^ 0.5) .* (eta.*(1 - eta)+1)  );
    epsilon_hat = (one_mat'.^3) ./ kappa ./ b .* (1 - b);
    q_hat = (2*k_hat).^(0.5);
    
    u_hat_ini = u_hat;
    c_hat_ini =  c_hat;
    c_hat_i_ini = c_hat_i;
    c_i_ini = c_i;
    c_ini = c;
    k_hat_ini = k_hat; 
    epsilon_hat_ini = epsilon_hat;
    q_hat_ini = q_hat;
    
    nu_t = C_mu .* k_hat .^ 2 ./ epsilon_hat;
    nu_tc = nu_t ./ sigma_sc; % was written as just "sigma"

    
    c_hat_i_int_ini = zeros(D_frac,1); % preallocating here is done by me!
    c_i_int_ini = zeros(D_frac,1);
    if (abs(sum(cb_i)) >= tol)
        for j = 1:D_frac
            [~, ~, ~, c_hat_i_int_ini(j)] = SimpsonInt( Nz, deta, c_hat_i_ini(:,j)', c_hat_i_int_ini(j) );
            [~, ~, ~, c_i_int_ini(j)] = SimpsonInt( Nz, deta, c_i_ini(:,j)', c_i_int_ini(j) );
        end   
    end
    c_hat_i_int = c_hat_i_int_ini;
    c_hat_int_ini = sum(c_hat_i_int_ini);
    c_i_int = c_i_int_ini;
    c_hat_int = c_hat_int_ini;
    c_int = sum(c_i_int);
    C_mean = c_int / (1-b);
    
    err = 10;
    i = 0;

    %% solution while loop
    while err > tol
        
        
        %% grab new values
        i = i + 1;
        u_hat_old = u_hat;
        c_hat_old = c_hat;
        c_hat_i_old = c_hat_i;
        c_old = c;
        c_i_old = c_i;
        k_hat_old = k_hat;
        epsilon_hat_old = epsilon_hat;
        q_hat_old = q_hat;
         
        
        %% calculate the new master length scale (DIMENSIONAL)
        var_1 = 0;
        var_2 = 0;
        [~,~,~,var_1] = SimpsonInt(Nz, dz, q_hat_old.*u_star, var_1);
        [~,~,~,var_2] = SimpsonInt(Nz, dz, q_hat_old.*u_star.*z, var_2);
        l_0 = alpha*var_2/var_1;
        l = l_0*kappa.*z ./ (kappa.*z+l_0);
        l_hat = l ./ h;
        
        
        %% calcualte the eddy viscosity profile
        dudz = (1-eta)./(nu_t+Re_flag*1/Re_tau);
        dcdz = zeros(Nz,1);
        dcdz_i = zeros(Nz,D_frac);
        for j = 1:D_frac
            % dcdz = dcdz + -V(j).*c_hat_i_old(:,j) ./ (nu_tc+Re_flag*1d0 / Re_tau / sigma_sc); % was "sigma", changed to "sigma_sc"
            dcdz_i(:, j) = -V(j).*c_hat_i_old(:,j) ./ (nu_tc+Re_flag*1d0 / Re_tau / sigma_sc);
        end
        dcdz = sum(dcdz_i,2);
        
        G_H = zeros(1,Nz);
        S_M = zeros(1,Nz);
        S_H = zeros(1,Nz);
        for j = 1:Nz
            if (q_hat_old(j) == 0d0)
                nu_t(j) = 0d0;
                nu_tc(j) = 0d0;
            else
                % Upper and lower bound limited
                G_H(j) = max( [Ri_flag*Ri_tau*((l_hat(j)/q_hat_old(j)).^2)*dcdz(j), -0.28] );
                G_H(j) = min( [G_H(j), 0.0223] ) ;
                % Upper and lower bound relaxed 
                S_M(j) = A_1*(1 - 3*C_1 - 6*A_1/B_1 - 3*A_2*G_H(j) * ((B_2 - 3*A_2) * (1 - 6*A_1/B_1) - 3*C_1*(B_2+6*A_1))) ...
                            /(1-3*A_2*G_H(j)*(6*A_1+B_2))/(1-9*A_1*A_2*G_H(j));
                S_H(j) = A_2*(1 - 6*A_1/B_1)/(1 - 3*A_2*G_H(j)*(6*A_1+B_2));
                nu_t(j) = q_hat_old(j)*l_hat(j)*S_M(j);
                nu_tc(j) = q_hat_old(j)*l_hat(j)*S_H(j);
            end
        end
        nu_t(1) = (1 - b - Re_flag / Re_tau*dudz(1)) / dudz(1);
        
        G_H_i = zeros(D_frac,Nz);
        S_M_i = zeros(D_frac,Nz);
        for d = 1:D_frac
            for j = 1:Nz
                if (q_hat_old(j) == 0d0)
                    nu_t_i(d,j) = 0d0;
                else
                    % Upper and lower bound limited
                    G_H_i(d,j) = max( [Ri_flag*Ri_tau*((l_hat(j)/q_hat_old(j)).^2)*dcdz_i(j,d), -0.28] );
                    G_H_i(d,j) = min( [G_H_i(d,j), 0.0223] ) ;
                    % Upper and lower bound relaxed 
                    S_M_i(d,j) = A_1*(1 - 3*C_1 - 6*A_1/B_1 - 3*A_2*G_H_i(d,j) * ((B_2 - 3*A_2) * (1 - 6*A_1/B_1) - 3*C_1*(B_2+6*A_1))) ...
                                /(1-3*A_2*G_H_i(d,j)*(6*A_1+B_2))/(1-9*A_1*A_2*G_H_i(d,j));
                    
                    nu_t_i(d,j) = q_hat_old(j)*l_hat(j)*S_M_i(d,j)*C_frac(d);
                end
            end
            nu_t_i(d, 1) = (1 - b - Re_flag / Re_tau*dudz(1)) / dudz(1);
        end
        
        
        %% diffusion for c
        for j = 1:D_frac
            % for the near bed specification case
            % for each grain size class
            if ((abs(V(j))>= tol) && (abs(C_bed_frac(j)) > 0d0))
                % if the settling velocity is faster than tol and frac>0
                SS = zero_mat;
                if (c_hat_i_old(Nz,j) <  1d-16) 
                    % if the surface conc is very small, do not worry about mass conservation
                    [c_hat_i_out] = UnsteadyAdvDiff1D_TVD(scheme, beta_TVD, Nz, eta, deta*one_mat, dt, f, -V(j)*one_mat, ... % missing deta*one_mat?
                                    nu_tc+Re_flag*1/Re_tau/sigma_sc, SS, c_hat_i_old(:, j)', c_hat_i(:, j)', ...
                                    1, -V(j), 0, 0, 1, c_hat_i_old(Nz,j), ...
                                    1, -V(j), 0, 0, 1, 0);
                    c_hat_i(:,j) = c_hat_i_out';
                    p1 = 0;
                    [~, ~, ~, p1] = SimpsonInt(Nz, deta, c_hat_i(:,j), p1);
                    c_hat_i(:,j) = c_hat_i(:,j) .* c_hat_i_int(j) ./ p1;
                else
                    % use a mass conservative boundary condition
                    [c_hat_i_out] = UnsteadyAdvDiff1D_consC(Nz, eta, deta, dt, f, -V(j)*one_mat, ...
                                            nu_tc'+Re_flag*1/Re_tau/sigma_sc, SS, ...
                                            c_hat_i_old(:,j)', c_hat_i(:,j)', c_hat_i_int(j), ...
                                            1, -V(j), 0, 0, 1, c_hat_i_old(Nz,j), ...
                                            1, -V(j), 0, 0, 1, c_hat_i_old(Nz,j));
                    c_hat_i(:,j) = c_hat_i_out';
                    if (c_hat_i(Nz,j) < 0d0)
                        c_hat_i(Nz,j) = 0d0;
                    end
                end
                % do a renormalization to ensure the mass conservation
                c_hat_i(:,j) = c_hat_i(:,j).*C_bed_frac(j).*sum(cb_i)./(c_hat_i(1,j)*C_mean);
            else
                % if frac ~ 0
                c_hat_i(:,j) = c_hat_i_old(:,j);
            end
        end
        % sum it up into a single c
        c_hat = sum(c_hat_i, 2);
        
        
        %% normalize c --- only if there is sediment
        if ~(C_mean==0)
            [~, ~, ~, c_hat_int] =  SimpsonInt(Nz, deta, c_hat, c_hat_int);
            C_mean_new = C_mean * c_hat_int / (1d0 - b);	% update C_mean
            Ri_tau = R * C_mean_new / S_0;		% update Ri_tau
            c_hat_i = (c_hat_i .* C_mean) ./ C_mean_new	;% renormalize c_i_new
            c_hat = (c_hat .* C_mean) ./ C_mean_new;		% renormalize c_new
            C_mean = C_mean_new;
            for j = 1:D_frac
                 [~, ~, ~, c_hat_i_int(j)] =  SimpsonInt(Nz, deta, c_hat_i(:,j), c_hat_i_int(j));
            end
            [~, ~, ~, c_hat_int] = SimpsonInt(Nz, deta, c_hat, c_hat_int);
            [~, ~, ~, c_int] =  SimpsonInt(Nz, dz, c, c_int);
        end        
        
        
        %% diffusion for u
        SS = one_mat;
        [u_hat] = SteadyAdvDiff1D_TVD(scheme,beta_TVD,Nz,eta,deta*one_mat,zero_mat,nu_t'+Re_flag*1./Re_tau,SS,u_hat', ...
                                          0d0,1d0,u_hat_old(1), ... %1/kappa*(log(b/k_s)+kappa*(8.5d0+(1d0/kappa*log(Re_tau*1*k_s)-3d0)*exp(-0.121d0*(log(Re_tau*1*k_s)).^2.42d0))), ...
                                          1d0,0d0,0d0);
        u_hat = u_hat';
        
        
        %% diffusion for tke balance
        if c_hat_i(1,:) > 0
            B_b = Ri_flag * Ri_tau * sum(V.*c_hat_i(1,:)/ sum(c_hat_i(1,:))); % some boundary thing?
        else
            B_b = (1.^2) ./ (C_mu .^ 0.5) .* (1 - b);
        end
        p1 = 1 / kappa / b-MO_a * B_b;

        SS = (2d0*nu_t.*dudz.^2) + (2d0*Ri_flag*Ri_tau*nu_tc.*dcdz) - (2*q_hat_old.^3 ./ B_1 ./ l_hat);
        [k_hat] = UnsteadyAdvDiff1D_TVD(0,     beta_TVD, Nz, z, dz*one_mat', dt, f, zero_mat', ...
                                        0.2*l_hat'.*q_hat_old'+Re_flag*1/Re_tau, SS', k_hat_old', k_hat', ...
                                        0d0, 1d0, k_hat_old(1), 1, 0, 0, ...
                                        0d0, 1d0, max([(B_1*(1-b-Re_flag/Re_tau*p1).*l_hat(1)*p1), 0]).^(2/3), 1, 0, 0);
        k_hat = k_hat';
        k_hat(k_hat < 0) = k_hat_old(k_hat < 0);
        q_hat = real(sqrt(k_hat));
        
        
        %% calculate an error and repeat
        err = max( [max( abs( (u_hat(err_idx)-u_hat_old(err_idx)) ./ (u_hat_old(err_idx)+1d-4) ) ), ...
                    max( abs( (c_hat(err_idx)-c_hat_old(err_idx)) ./ (c_hat_old(err_idx)+1d-4) ) ), ...
                    max( abs( (q_hat(err_idx)-q_hat_old(err_idx)) ./ (q_hat_old(err_idx)+1d-4) ) )] );
        
        if i > maxiter
            err = 0;
        end
          
        
        %% print update and plot
        if (mod(i, parkat) == 0)
            if verbose
                disp( ['                iteration: ', num2str(i), ', near-bed conc: ', num2str(c_hat(1)), ...
                       ', depth avg conc: ', num2str(c_int), ', depth avg vel: ', num2str(mean(u_hat).*u_star), ', error: ', num2str(err)] )
            end
            
            % plot
            if solution_plot
                subplot(1,4,1); hold on;
                    cla
                    plot(u_non, eta, 'k--')
                    plot(u_hat_old, eta, 'k')
                    plot(u_hat, eta, 'r')
                    xlabel('u/u_*')
                    ylabel('z/h')
                subplot(1,4,2); hold on;
                    cla
                    plot(c_non.*sum(cb_i), eta, 'k--')
                    plot(c_hat.*C_mean, eta, 'r')
                    xlabel('c')
                subplot(1,4,3); hold on
                    cla
                    plot(nu_t_non, eta, 'k--')
                    plot(nu_t, eta, 'r')
                    xlabel('\nu_t')
                subplot(1,4,4); hold on;
                    cla
                    plot(k_hat_ini, eta, 'k:')
                    plot(k_hat, eta, 'r')
                    xlabel('k/u_*^2')
                drawnow
            end
        end

    end % end while loop for main solution
    
    
    %% finalize and make the plot
    if final_plot
        % plot
        close all
        figure()
        subplot(1,4,1); hold on;
            cla
            plot(simple.u_non, eta, 'k--')
            plot(u_hat_old, eta, 'k')
            plot(u_hat, eta, 'r')
            xlabel('u/u_*')
            ylabel('z/h')
        subplot(1,4,2); hold on;
            cla
            plot(simple.c_non.*sum(cb_i), eta, 'k--')
            plot(c_hat.*C_mean, eta, 'r')
            xlabel('c')
        subplot(1,4,3); hold on
            cla
            plot(simple.nu_t_non, eta, 'k--')
            plot(nu_t, eta, 'r')
            xlabel('\nu_t')
        subplot(1,4,4); hold on;
            cla
            plot(k_hat_ini, eta, 'k:')
            plot(k_hat, eta, 'r')
            xlabel('k/u_*^2')
        drawnow
    end
    
    %% pack object to return out
    strat.u_hat = u_hat;
    strat.Us = u_hat .* u_star;
    strat.c_hat = c_hat;
    strat.CsSum = c_hat .* C_mean .* 2650;
    strat.k_hat = k_hat;
    strat.k = k_hat .* (u_star * u_star);
    strat.l_0 = l_hat .* h;
    strat.Cs = c_hat_i .* C_mean .* 2650;
    strat.Zs = z;
    strat.ZsSum = z;
    strat.eta = eta;
    strat.C_mean = C_mean;
    strat.nu_t = nu_t;
    strat.nu_tc = nu_tc;
    strat.nu_t_i = nu_t_i;
    
end