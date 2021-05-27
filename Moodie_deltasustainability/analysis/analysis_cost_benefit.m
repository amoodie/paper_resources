function analysis_cost_benefit()

    clear all
    
    %% load data
    % add things to path
    source_path = genpath(fullfile('..', 'source'));
    input_path = genpath(fullfile('..', 'input_data'));
    output_data_path = genpath(fullfile('..', 'output_data'));
    home_path = genpath(fullfile('~', 'Documents', 'MATLAB'));
    addpath(source_path, input_path, output_data_path, plotting_path, home_path);
    [CP] = color_palette();
    label_interpreter = 'latex';
    
    printOut = true;
    print_root = '.';
    bgcolor = [1 1 1]; % [0.94 0.94 0.94];
    mksize = 40;
    
    % load engineered data
    load(fullfile('engineered', 'engineered_single_avulsion_location_list.mat'), 'avulsion_location_list');
    load(fullfile('engineered_single_analysis_summary.mat')) % loads eng into mem
    
    
    % load engineered spinup data
    spin = load(fullfile('engineered', 'spinups', ['engineered_spinup_data_single.mat']));
    spin.analysis = load(fullfile('engineered', 'spinups', ['engineered_spinup_analysis_single.mat']));
    spin.analysis = spin.analysis.analysis;
    spin.avul_time = spin.analysis.avul_time_mean;
    spin.avul_len = spin.analysis.avul_len_mean;
    spin.lobe_len = spin.analysis.lobe_len_mean;
    

    %% cost functions
    %costfn_donothing_nonannualized = @(beta, theta, lmb_f, lmb_l, R_l) ((pi/4) + 0.5.*beta.*R_l.*lmb_l - tand(theta/2).*(1-beta).^2.*(1+(1./lmb_f)));
    %costfn_donothing = @(beta, theta, lmb_f, lmb_l, R_l) ((pi/4) + 0.5.*beta.*R_l.*lmb_l - tand(theta/2).*(1-beta).^2.*(1+(1./lmb_f)));
    %costfn_doavulsion_nonannualized = @(L_A_star, L_f_star, L_d_star, T_A_star, alpha, beta, theta, R_l, lmb_f, lmb_l, lmb_A) (T_A_star.*((pi/4).*L_d_star.^2 + 0.5.*beta.*R_l.*lmb_l) - lmb_A.*(1+alpha.*L_f_star) - (tand(theta/2).*(L_f_star).^2).*(1+(1./lmb_f)));
    %costfn_doavulsion = @(L_A_star, L_f_star, L_d_star, T_A_star, alpha, beta, theta, R_l, lmb_f, lmb_l, lmb_A) (T_A_star.*((pi/4).*L_d_star.^2 + 0.5.*beta.*R_l.*lmb_l) - lmb_A.*(1+alpha.*L_f_star) - (tand(theta/2).*(L_f_star).^2).*(1+(1./lmb_f)))./T_A_star;

    costfn_donothing = @(beta, theta, lmb_f, lmb_l, R_l) ...
                         1 + (0.5.*beta.*R_l.*lmb_l.*(1./((pi/4)))) - (((tand(theta/2).*(1-beta).^2)).*(1+lmb_f));
    costfn_doavulsion = @(L_D_star, L_f_star, L_d_star, T_A_star, ...
                          alpha, beta, theta, R_l, lmb_f, lmb_l, lmb_D) ...
                          1 - ((lmb_D./T_A_star).*(1+alpha.*L_f_star)) + ...
                          (0.5.*beta.*R_l.*lmb_l.*(1./((pi/4).*L_d_star.^2))) - ...
                          (((tand(theta/2).*(L_f_star).^2)./T_A_star).*(1+lmb_f));
                      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These cost functions below assume that flooding area is a rectangle 0.4 Lb in width.
%     costfn_donothing = @(beta, theta, lmb_f, lmb_l, R_l) ...
%                         1 + (0.5.*beta.*R_l.*lmb_l.*(1./((pi/4)))) - (((0.05.*(1-beta))).*(1+lmb_f));
%     costfn_doavulsion = @(L_A_star, L_f_star, L_d_star, T_A_star, ...
%                          alpha, beta, theta, R_l, lmb_f, lmb_l, lmb_A) ...
%                          1 - ((lmb_A./T_A_star).*(1+alpha.*L_f_star)) + ...
%                          (0.5.*beta.*R_l.*lmb_l.*(1./((pi/4).*L_d_star.^2))) - ...
%                          (((0.05.*(L_f_star))./T_A_star).*(1+lmb_f));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % split doavulsion into each term
    costfn_doavulsion_lobe = @(L_D_star, L_f_star, L_d_star, T_A_star, ...
                          alpha, beta, theta, R_l, lmb_f, lmb_l, lmb_D) ...
                          (0.5.*beta.*R_l.*lmb_l.*(1./((pi/4).*L_d_star.^2)));
    costfn_doavulsion_avul = @(L_D_star, L_f_star, L_d_star, T_A_star, ...
                          alpha, beta, theta, R_l, lmb_f, lmb_l, lmb_D) ...
                          - ((lmb_D./T_A_star).*(1+alpha.*L_f_star));
    costfn_doavulsion_flood = @(L_D_star, L_f_star, L_d_star, T_A_star, ...
                          alpha, beta, theta, R_l, lmb_f, lmb_l, lmb_D) ...
                          - (((tand(theta/2).*(L_f_star).^2)./T_A_star).*(1+lmb_f));

    %% avulsion timescale lengthscale relationship functions 
    % L_A_star to T_A_star functions
    return_fn_squared = @(L_D_star) L_D_star.^2;
    return_fn_squareroot = @(L_D_star) L_D_star.^(0.5);
    return_fn_linear = @(L_D_star) L_D_star;
    
    
    %% hydraulic assumptions
    natural.L_D_star = 1.0; % avulsion length scale if we do nothing
    natural.T_A_star = 1.0; % avulsion time scale if we do nothing
    natural.lambda_A = 0.0; % avulsion cost if we do nothing
    
    
    %% define montecarlo distributions and params
    
    % distribution parameters
    lmb_l_ms = [1, 0.1]; % normal, mean std
    lmb_f_ms = [2, 0.1]; % normal, mean std
    lmb_D_ms = [0.01, 0.03]; % normal, mean std
    beta_ms = [0.4, 0.6]; % uniform, limits
    alpha_ms = [0.01, 0.05]; % uniform, limits
    R_l_ms = [0.22, 0.08]; % normal, mean std
    % theta_ms = [20, 8]; % normal, mean std
    theta_ms = [2, 6]; % gamma, shape a, b
    
    % make distributions
    lmb_l_dist = truncate(makedist('Normal', lmb_l_ms(1), lmb_l_ms(2)), 0, inf);
    lmb_f_dist = truncate(makedist('Normal', lmb_f_ms(1), lmb_f_ms(2)), 0, inf);
    lmb_D_dist = truncate(makedist('Normal', lmb_D_ms(1), lmb_D_ms(2)), 0, inf);
    beta_dist = makedist('Uniform', beta_ms(1), beta_ms(2));
    alpha_dist = makedist('Uniform', alpha_ms(1), alpha_ms(2));
    R_l_dist = truncate(makedist('Normal', R_l_ms(1), R_l_ms(2)), 0, inf);
    theta_dist = truncate(makedist('Gamma', theta_ms(1), theta_ms(2)), 0, inf);
    
    % sampling functions
    sample_lmb_l = @(n) random(lmb_l_dist,n,1);        % @() lmb_l_ms(1) + lmb_l_ms(2).*randn(1);
    sample_lmb_f = @(n) random(lmb_f_dist,n,1);        % @() lmb_f_ms(1) + lmb_f_ms(2).*randn(1);
    sample_lmb_D = @(n) random(lmb_D_dist,n,1);        % @() lmb_A_ms(1) + lmb_A_ms(2).*randn(1);
    sample_beta =  @(n) random(beta_dist,n,1);         % @() beta_ms(1) + (beta_ms(2)-beta_ms(1))*rand(n,1);
    sample_alpha = @(n) random(alpha_dist,n,1);        % @() alpha_ms(1) + (alpha_ms(2)-alpha_ms(1))*rand(n,1);
    sample_R_l =   @(n) random(R_l_dist,n,1);          % @() R_l_ms(1) + R_l_ms(2).*randn(1);
    sample_theta = @(n) random(theta_dist,n,1);        % @() theta_ms(1) + theta_ms(2).*randn(1);
    
    L_D_star_list = linspace(0.01, 1.8, 50);
    n_mc = 1000;
    lamb_fs_plot = [1, 5, 10]; % explicit ones to plot and eval over
    lamb_Ds_plot = [0.01, 0.1, 20]; % explicit ones to plot and eval over
    
    %% montecarlo for different avulsion TA~LA functions
    mced_func_donothing = NaN(n_mc,length(L_D_star_list),3);
    mced_func_doavulsion = NaN(n_mc,length(L_D_star_list),3);
    mced_func_sustnum = NaN(n_mc,length(L_D_star_list),3);
    mced_func_doavulsion_parts = NaN(n_mc,length(L_D_star_list),3,3);
    
    smpl_lmb_l = sample_lmb_l(n_mc);
    smpl_lmb_f = sample_lmb_f(n_mc);
    smpl_lmb_D = sample_lmb_D(n_mc);
    smpl_beta = sample_beta(n_mc);
    smpl_R_l = sample_R_l(n_mc);
    smpl_alpha = sample_alpha(n_mc);
    smpl_theta = sample_theta(n_mc);
    for j = 1:length(L_D_star_list)
        L_D_star_j = L_D_star_list(j);
        % set the model (evaluated into cols, as sqrd, sqrt, lin
        T_A_star_j_sqrd = return_fn_squared(L_D_star_j);
        T_A_star_j_sqrt = return_fn_squareroot(L_D_star_j);
        T_A_star_j_linear = return_fn_linear(L_D_star_j);
        for k = 1:n_mc
            
            L_d_star_ijk = 1; % max(1, L_A_star_j);
            L_f_star_ijk = return_L_f_star(L_D_star_j, smpl_beta(k), smpl_R_l(k));

            mced_func_donothing(k,j,1:3) = costfn_donothing(smpl_beta(k), smpl_theta(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_R_l(k));
            mced_func_doavulsion(k,j,1) = costfn_doavulsion(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_sqrd, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_func_doavulsion(k,j,2) = costfn_doavulsion(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_sqrt, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_func_doavulsion(k,j,3) = costfn_doavulsion(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_linear, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            
            mced_func_sustnum(k,j,1:3) = mced_func_doavulsion(k,j,1:3) ./ mced_func_donothing(k,j,1:3);
            
            mced_func_doavulsion_parts(k,j,1,1) = costfn_doavulsion_lobe(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_sqrd, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_func_doavulsion_parts(k,j,2,1) = costfn_doavulsion_lobe(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_sqrt, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_func_doavulsion_parts(k,j,3,1) = costfn_doavulsion_lobe(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_linear, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_func_doavulsion_parts(k,j,1,2) = costfn_doavulsion_avul(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_sqrd, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_func_doavulsion_parts(k,j,2,2) = costfn_doavulsion_avul(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_sqrt, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_func_doavulsion_parts(k,j,3,2) = costfn_doavulsion_avul(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_linear, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_func_doavulsion_parts(k,j,1,3) = costfn_doavulsion_flood(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_sqrd, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_func_doavulsion_parts(k,j,2,3) = costfn_doavulsion_flood(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_sqrt, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_func_doavulsion_parts(k,j,3,3) = costfn_doavulsion_flood(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j_linear, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
        end
    end
    
    %% evaluate the cost for yellow river
    % reuse the random samples from above
    mced_YR_donothing = NaN(n_mc,length(eng.tbl.avulloc));
    mced_YR_doavulsion = NaN(n_mc,length(eng.tbl.avulloc));
    mced_YR_sustnum = NaN(n_mc,length(eng.tbl.avulloc));
    for j = 1:size(eng.tbl,1)
        L_D_star_j = eng.tbl.avulloc(j);
        for k = 1:n_mc
            good_YR_T_A_star = false;
            cnt = 0;
            while ~good_YR_T_A_star
                YR_T_A_star = (eng.tbl.TA_mean(j) + eng.tbl.TA_std(j).*randn(1)) ./ spin.avul_time; % do a random draw and then nromalize it for YR T_A_star
                %YR_T_A_star = max(0, YR_T_A_star);
                if YR_T_A_star > 0
                    good_YR_T_A_star = true;
                else
                    cnt = cnt + 1;
                end
            end
            L_d_star_ijk = 1; % max(1, L_A_star_j);
            L_f_star_ijk = return_L_f_star(L_D_star_j, smpl_beta(k), smpl_R_l(k));

            mced_YR_donothing(k,j) = costfn_donothing(smpl_beta(k), smpl_theta(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_R_l(k));
            mced_YR_doavulsion(k,j) = costfn_doavulsion(L_D_star_j, L_f_star_ijk, L_d_star_ijk, YR_T_A_star, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), smpl_lmb_f(k), smpl_lmb_l(k), smpl_lmb_D(k));
            mced_YR_sustnum(k,j) = mced_YR_doavulsion(k,j) ./ mced_YR_donothing(k,j);
        end
    end
    
    
    
    %% montecarlo for sqrt function only varying lambda f
    % reuse the random samples from above
    mced_lmb_f_donothing = NaN(n_mc,length(L_D_star_list),length(lamb_fs_plot));
    mced_lmb_f_doavulsion = NaN(n_mc,length(L_D_star_list),length(lamb_fs_plot));
    mced_lmb_f_sustnum = NaN(n_mc,length(L_D_star_list),length(lamb_fs_plot));
    for i = 1:length(lamb_fs_plot)
        lamb_fs_i = lamb_fs_plot(i);
        for j = 1:length(L_D_star_list)
            L_D_star_j = L_D_star_list(j);
            T_A_star_j = return_fn_squareroot(L_D_star_j);
            for k = 1:n_mc
                L_d_star_ijk = 1;
                L_f_star_ijk = return_L_f_star(L_D_star_j, smpl_beta(k), smpl_R_l(k));
                
                mced_lmb_f_donothing(k,j,i) = costfn_donothing(smpl_beta(k), smpl_theta(k), lamb_fs_i, smpl_lmb_l(k), smpl_R_l(k));
                mced_lmb_f_doavulsion(k,j,i) = costfn_doavulsion(L_D_star_j, L_f_star_ijk, L_d_star_ijk, T_A_star_j, smpl_alpha(k), smpl_beta(k), smpl_theta(k), smpl_R_l(k), lamb_fs_i, smpl_lmb_l(k), smpl_lmb_D(k));
                mced_lmb_f_sustnum(k,j,i) = mced_lmb_f_doavulsion(k,j,i) ./ mced_lmb_f_donothing(k,j,i);
            end
        end
    end
    
    %% MAIN figure init
    MAINfig = figure();
    set(MAINfig, 'Pos', [700 500 634 645]);

    %% natural avulsion profit parameter space
    % evaluation of natural profit for various lambda_f
    nat_lnes = NaN(length(lamb_fs_plot), 1);
    for i = 1:length(lamb_fs_plot)
        nat_lnes(i) = costfn_donothing(mean(beta_dist), mean(theta_dist), ...
                                       lamb_fs_plot(i), mean(lmb_l_dist), mean(R_l_dist));
    end
   
    % parameter space grid
    lin_lmb_l = linspace(0, 4, 50);
    lin_lmb_f = linspace(0, 10, 50);
    [mg_lmb_f, mg_lmb_l] = meshgrid(lin_lmb_f, lin_lmb_l);
    
    nat_mg = costfn_donothing(mean(beta_dist), mean(theta_dist), mg_lmb_f, mg_lmb_l, mean(R_l_dist));
    nat_mg = nat_mg;
    
    subplot(2,2,1); hold on;
        plot([0, 10], [1, 1], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        plot([1, 1], [0, 4], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        [C,h] = contour(mg_lmb_f, mg_lmb_l, nat_mg, 'LineWidth', 2);
        view([0 90])
        clabel(C,h, 'interpreter', label_interpreter, 'Color', [0.5 0.5 0.5])
        %colormap(ndvi_colormap)
        caxis([0.8 1.2])
    text(0.85, 0.85, 'a', 'units', 'normalized', 'FontSize', 12, 'BackgroundColor', [1 1 1])
    text(0.35, 0.92, '(natural)', 'units', 'normalized', 'FontSize', 10, 'BackgroundColor', [1 1 1])
    xticks([0:2:10])
    xlabel('flooding cost ($\lambda_f$)', 'interpreter', label_interpreter)
    ylabel('lobe land value ($\lambda_l$)', 'interpreter', label_interpreter)
    axis square
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
    box on
    set(gca, 'Layer', 'top')
    
    
    %% artificial avulsion profit parameter space
    % parameter space grid
    lin_lmb_D = logspace(-4, -1, 50);
    [mg_lmb_f, mg_lmb_D] = meshgrid(lin_lmb_f, lin_lmb_D);
    
    L_D_test = 1;
    art_mg = costfn_doavulsion(L_D_test, return_L_f_star(L_D_test, mean(beta_dist), mean(R_l_dist)), L_d_star_ijk, ...
                               return_fn_squareroot(L_D_test), mean(alpha_dist), mean(beta_dist), mean(theta_dist), ...
                               mean(R_l_dist), mg_lmb_f, mean(lmb_l_dist), mg_lmb_D);
    art_mg = art_mg;
    
    subplot(2,2,2); hold on;
        plot([1, 1], [1e-4, 1e-1], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        [C,h] = contour(mg_lmb_f, mg_lmb_D, art_mg, 'LineWidth', 2);
        view([0 90])
        clabel(C,h, 'interpreter', label_interpreter, 'Color', [0.5 0.5 0.5])
        set(gca, 'yscale', 'log')
        %colormap(ndvi_colormap)
        caxis([0.8 1.2])
    text(0.85, 0.85, 'b', 'units', 'normalized', 'FontSize', 12, 'BackgroundColor', [1 1 1])
    text(0.35, 0.92, '(artificial)', 'units', 'normalized', 'FontSize', 10, 'BackgroundColor', [1 1 1])
    xticks([0:2:10])
    xlabel('flooding cost ($\lambda_f$)', 'interpreter', label_interpreter)
    ylabel('diversion cost ($\lambda_D$)', 'interpreter', label_interpreter)
    axis square
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
    box on
    set(gca, 'Layer', 'top')
   

    
    %% plot three function models mced
    pcmap = summer(size(eng.tbl,1));
    pcmap = repmat(pcmap(12,:), size(eng.tbl,1), 1);
    mced_func_mean = squeeze(mean(mced_func_doavulsion, 1));
    mced_func_std = squeeze(std(mced_func_doavulsion, [], 1));
    mced_func_95 = (1.96*mced_func_std);
    mced_YR_mean = squeeze(nanmean(mced_YR_doavulsion, 1));
    mced_YR_std = squeeze(std(mced_YR_doavulsion, [], 1));
    mced_YR_sustnmum_mean = squeeze(nanmean(mced_YR_sustnum, 1));
    mced_YR_sustnmum_std = squeeze(std(mced_YR_sustnum, [], 1));
    mced_func_sustnum_mean = squeeze(nanmean(mced_func_sustnum, 1));
    mced_func_sustnum_std = squeeze(std(mced_func_sustnum, [], 1));
    mced_func_sustnum_95 = (1.96*mced_func_sustnum_std);
    [mced_func_mean_maxv, mced_func_mean_maxi] = max(mced_func_mean,[],1);
    [mced_func_sustnum_mean_maxv, mced_func_sustnum_mean_maxi] = max(mced_func_sustnum_mean,[],1);
    
    func_cmap.fills = [CP.c3_s2; CP.c2_s2; CP.c1_s2];
    func_cmap.lines = [CP.c3_s1; CP.c2_s1; CP.c1_s1];
    sustnum_ylims = [0.5, 1.4];
    
    figure(MAINfig);
    subplot(2, 2, 3); hold on
    plot([0, 2], [0, 0], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
    plot_as_sustnum = true;
    if ~plot_as_sustnum
        for i = 3:3
            fill([L_D_star_list, fliplr(L_D_star_list)], [mced_func_mean(:,i)+mced_func_95(:,i); flipud(mced_func_mean(:,i)-mced_func_95(:,i))], ...
                func_cmap.fills(i,:), 'EdgeColor', factor_color(func_cmap.lines(i,:), 1.1), 'FaceAlpha', 1)
        end
        for i = 3:3
            plot(L_D_star_list, mced_func_mean(:,i), 'Color', func_cmap.lines(i,:), 'LineWidth', 1.2)
        end
        errorbar(eng.tbl.avulloc, mced_YR_mean, mced_YR_std, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)   
        scatter(eng.tbl.avulloc, mced_YR_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
        ylim([0.5, 1.25])        
        yticks(0.5:0.25:1.25)
        ylabel('${\lambda_\Pi}/{T_A}^*$')
    else  % plot as sustnum = true
        for i = 2:2
            fill([L_D_star_list, fliplr(L_D_star_list)], [mced_func_sustnum_mean(:,i)+mced_func_sustnum_95(:,i); flipud(mced_func_sustnum_mean(:,i)-mced_func_sustnum_95(:,i))], ...
                func_cmap.fills(i,:), 'EdgeColor', factor_color(func_cmap.lines(i,:), 1.1), 'FaceAlpha', 0.7)
        end
        for i = 2:2
            plot(L_D_star_list, mced_func_sustnum_mean(:,i), 'Color', func_cmap.lines(i,:), 'LineWidth', 1.2)
            plot(L_D_star_list(mced_func_sustnum_mean_maxi(i)), mced_func_sustnum_mean_maxv(i), 'o', ...
                'MarkerFaceColor', factor_color(func_cmap.lines(i,:), 1), 'MarkerEdgeColor', 'none', 'MarkerSize', 7)
        end
        plot([0, 2], [1, 1], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        errorbar(eng.tbl.avulloc, mced_YR_sustnmum_mean, mced_YR_sustnmum_std, 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 0)   
        scatter(eng.tbl.avulloc, mced_YR_sustnmum_mean, mksize, pcmap, 'filled', 'MarkerEdgeColor', 'k')
        ylim(sustnum_ylims)
        %yticks([0.5, 0.75, 1, 1.15, 1.3])
        ylabel('sustainability number ($\mathcal{S}$)')
    end
    xlabel('diversion length (${L_D}^*$)')
    text(0.85, 0.85, 'c', 'units', 'normalized', 'FontSize', 12)
    xlim([0 1.7])
    axis square
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
    box on
    set(gca, 'Layer', 'top')
    
    
    
    %% plot sqrt model for different lambda_f values
    mced_lmb_f_mean = squeeze(mean(mced_lmb_f_doavulsion, 1));
    mced_lmb_f_std = squeeze(std(mced_lmb_f_doavulsion, [], 1));
    mced_lmb_f_dn_mean = squeeze(mean(mced_lmb_f_donothing, 1));
    mced_lmb_f_dn_std = squeeze(std(mced_lmb_f_donothing, [], 1));
    [mced_lmb_f_mean_maxv, mced_lmb_f_mean_maxi] = max(mced_lmb_f_mean./nat_lnes',[],1);
    
    lmbf_cmap = parula(length(lamb_fs_plot));
    figure(MAINfig);
    subplot(2, 2, 4); hold on
    plot_as_sustnum = true;
    if ~plot_as_sustnum
        for i = 1:length(lamb_fs_plot)
            f = fill([L_D_star_list, fliplr(L_D_star_list)], [mced_lmb_f_mean(:,i)+mced_lmb_f_std(:,i); flipud(mced_lmb_f_mean(:,i)-mced_lmb_f_std(:,i))], ...
                lmbf_cmap(i,:), 'EdgeColor', factor_color(lmbf_cmap(i,:), 1.1), 'FaceAlpha', 0.2);
        end
        for i = 1:length(lamb_fs_plot)
            plot(L_D_star_list, mced_lmb_f_mean(:,i), '-', 'Color', factor_color(lmbf_cmap(i,:),0.8), 'LineWidth', 1.2)
            text(L_D_star_list(11), mced_lmb_f_mean(11,i), ['$\lambda_f = ', num2str(lamb_fs_plot(i)), '$'], 'Color', lmbf_cmap(i,:))
            plot(repmat([0; 2], 1, length(nat_lnes)), [nat_lnes, nat_lnes]', ':', ...
                'Color', factor_color(lmbf_cmap(i,:), 0.8), 'LineWidth', 1.2)
        end
        ylim([0.5, 1.25])
        yticks(0.5:0.25:1.25)
        ylabel('${\lambda_\Pi}/{T_A}^*$')
    else
        lamb_fs_label_pos = {[0.26 0.89, 10], [0.2 1.12, 12], [0.22 1.20 15]}; % [x, y, rotation]
        for i = 1:length(lamb_fs_plot)
            fill([L_D_star_list, fliplr(L_D_star_list)], [mced_lmb_f_mean(:,i)+mced_lmb_f_std(:,i); flipud(mced_lmb_f_mean(:,i)-mced_lmb_f_std(:,i))]./nat_lnes(i), ...
                lmbf_cmap(i,:), 'EdgeColor', factor_color(lmbf_cmap(i,:), 1.1), 'FaceAlpha', 0.3);
        end
        for i = 1:length(lamb_fs_plot)
            plot(L_D_star_list, mced_lmb_f_mean(:,i)./nat_lnes(i), '-', 'Color', factor_color(lmbf_cmap(i,:),0.8), 'LineWidth', 1.2)
            text(lamb_fs_label_pos{i}(1), lamb_fs_label_pos{i}(2), ['$\lambda_f = ', num2str(lamb_fs_plot(i)), '$'], ...
                'Color', factor_color(lmbf_cmap(i,:), 0.8), 'rotation', lamb_fs_label_pos{i}(3))
            plot(L_D_star_list(mced_lmb_f_mean_maxi(i)), mced_lmb_f_mean_maxv(i), 'o', ...
                'MarkerFaceColor', factor_color(lmbf_cmap(i,:),0.8), 'MarkerEdgeColor', 'none', 'MarkerSize', 7)
        end
        plot([0, 2], [1, 1], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        ylim(sustnum_ylims)
        ylabel('sustainability number ($\mathcal{S}$)')
    end
    xlabel('diversion length (${L_D}^*$)')
    text(0.85, 0.85, 'd', 'units', 'normalized', 'FontSize', 12)
    xlim([0 1.7])
    axis square
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
    box on
    set(gca, 'Layer', 'top')
    
    set(MAINfig, 'Pos', [700 500 570 552]);
    set(MAINfig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    if printOut
        print('-dpng', '-r300', [print_root, 'engineered_economics_plots.png']);
        set(gcf, 'InvertHardcopy', 'on')
        print('-depsc', '-opengl', '-r300', [print_root, 'engineered_economics_plots.eps']);
        pause(1);
        export_fig([print_root, 'engineered_economics_plots.pdf'])
    end
    

    %% plot three models and three mced evaluations
    fig = figure();
    subplot(1, 2, 1)
        hold on
        plot([0 1.7], [1, 1], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        plot((0:0.05:1.9), (0:0.05:1.9), '-', 'Color', CP.c1_s3, 'LineWidth', 3)
        text(0.46, 0.30, '${L_D}^*$', 'rotation', 55, 'Color', CP.c1_s1, 'fontweight', 'bold')
        plot((0:0.05:1.9), (0:0.05:1.9).^0.5, '-', 'Color', CP.c2_s3, 'LineWidth', 3)
        text(0.25, 0.68, '$\sqrt{{L_D}^*}$', 'rotation', 32, 'Color', CP.c2_s1, 'fontweight', 'bold', 'fontsize', 12)
        plot((0:0.05:1.9), (0:0.05:1.9).^2, '-', 'Color', CP.c3_s3, 'LineWidth', 3)
        text(0.72, 0.28, '$({L_D}^*)^{2}$', 'rotation', 60, 'Color', CP.c3_s1, 'fontweight', 'bold')
        text(0.85, 0.1, 'a', 'units', 'normalized', 'FontSize', 12)
        ylabel('time to subsequent avulsion (${{T_A}^*}$)', 'interpreter', label_interpreter)
        yticks(0:0.2:1.6)
        xlim([0 1.7])
        ylim([0 1.6])
        box on
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        set(gca, 'Layer', 'top')
        axis square
    subplot(1, 2, 2); hold on
        for i = 1:3
            fill([L_D_star_list, fliplr(L_D_star_list)], [mced_func_sustnum_mean(:,i)+mced_func_sustnum_95(:,i); flipud(mced_func_sustnum_mean(:,i)-mced_func_sustnum_95(:,i))], ...
                func_cmap.fills(i,:), 'EdgeColor', factor_color(func_cmap.lines(i,:), 1.1), 'FaceAlpha', 0.4)
        end
        for i = 1:3
            plot(L_D_star_list, mced_func_sustnum_mean(:,i), 'Color', func_cmap.lines(i,:), 'LineWidth', 1.2)
            plot(L_D_star_list(mced_func_sustnum_mean_maxi(i)), mced_func_sustnum_mean_maxv(i), 'o', ...
                'MarkerFaceColor', factor_color(func_cmap.lines(i,:), 1.1), 'MarkerEdgeColor', 'none', 'MarkerSize', 7)
        end
        plot([0, 2], [1, 1], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        ylim(sustnum_ylims)
        ylabel('sustainability number ($\mathcal{S}$)')
        xl = xlabel('diversion length (${L_D}^* = L_D/L_b$)', 'interpreter', label_interpreter);
        text(0.85, 0.85, 'b', 'units', 'normalized', 'FontSize', 12)
        xlim([0 1.7])
        axis square
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
        xl.Position(1) = xl.Position(1) - 1.25;
        xl.Position(2) = xl.Position(2) - 0.075;
    set(fig, 'Pos', [200 500 570 300]);
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    if printOut
        print('-dpng', '-r300', [print_root, 'engineered_threerelations.png']);
        set(gcf, 'InvertHardcopy', 'on')
        print('-depsc', '-opengl', '-r300', [print_root, 'engineered_threerelations.eps']);
        pause(1);
    end
    
    
    %% parts of cost function plot
    mced_func_doavulsion_parts_mean = squeeze(mean(mced_func_doavulsion_parts, 1));
    mced_func_doavulsion_parts_std = squeeze(std(mced_func_doavulsion_parts, [], 1));
    
    cat_cmap = lines(3);
    label_list = {'a', 'b', 'c'};
    name_list = {'${T_A}^*={{L_D}^*}^2$', '${T_A}^*=\sqrt{{L_D}^*}$', '${T_A}^*={{L_D}^*}$'};
    fig = figure();
    for f = 1:3
        subplot(1, 3, f); hold on;
        for i = 1:3
            fill([L_D_star_list, fliplr(L_D_star_list)], [mced_func_doavulsion_parts_mean(:,f,i)+mced_func_doavulsion_parts_std(:,f,i); flipud(mced_func_doavulsion_parts_mean(:,f,i)-mced_func_doavulsion_parts_std(:,f,i))], ...
                cat_cmap(i,:), 'EdgeColor', factor_color(cat_cmap(i,:), 1.1), 'FaceAlpha', 0.2);
        end
        for i = 1:3
            lne(i) = plot(L_D_star_list, mced_func_doavulsion_parts_mean(:,f,i), 'LineWidth', 1.2, 'Color', cat_cmap(i,:));
        end
        sm = plot(L_D_star_list, sum(mced_func_doavulsion_parts_mean(:,f,:),3), 'k-', 'LineWidth', 1.2);
        plot([0, 1.8], [0, 0], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        ylim([-0.4, 0.2])
        yticks(-0.6:0.2:0.4)
        if f == 1
            ylabel('term of ${\lambda_\Pi}/{T_A}^*$')
        end
        xlabel('${L_D}^*$')
        text(0.90, 0.90, label_list{f}, 'units', 'normalized', 'FontSize', 12)
        text(0.10, 0.90, name_list{f}, 'units', 'normalized', 'FontSize', 10)
        xlim([0 1.7])
        axis square
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
    end
    legend([lne'; sm], {'lobe revenue', 'diversion cost', 'flooding cost', 'sum of terms'}, 'location', 'southeast', 'interpreter', 'latex')
    set(fig, 'Pos', [700 500 815 285]);
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    if printOut
        print('-dpng', '-r300', [print_root, 'engineered_doavulsion_parts.png']);
        set(gcf, 'InvertHardcopy', 'on')
        print('-depsc', '-opengl', '-r300', [print_root, 'engineered_doavulsion_parts.eps']);
        pause(1);
    end
    
     %% sustainability number parameter space
    % parameter space grids
    lin_lmb_D = logspace(-5, -1, 50);
    lin_lmb_f = linspace(0, 10, 50);
    lin_lmb_l = linspace(0, 10, 50);
    lin_alpha = linspace(0.01, 0.1, 50);
    lin_beta = linspace(0.4, 0.6, 50);
    lin_theta = linspace(1, 45, 50);
    lin_R_l = linspace(0.05, 0.5, 50);
    lin_L_D_star = linspace(0.05, 1.6, 50);
    [mg_L_D_star, mg_lmb_f] = meshgrid(lin_L_D_star, lin_lmb_f);
    [~, mg_lmb_D] = meshgrid(lin_L_D_star, lin_lmb_D);
    [~, mg_lmb_l] = meshgrid(lin_L_D_star, lin_lmb_l);
    [~, mg_alpha] = meshgrid(lin_L_D_star, lin_alpha);
    [~, mg_beta] = meshgrid(lin_L_D_star, lin_beta);
    [~, mg_theta] = meshgrid(lin_L_D_star, lin_theta);
    [~, mg_R_l] = meshgrid(lin_L_D_star, lin_R_l);
    
    sust_da_mg_lmb_D = costfn_doavulsion(mg_L_D_star, return_L_f_star(mg_L_D_star, mean(beta_dist), mean(R_l_dist)), L_d_star_ijk, ...
                               return_fn_squareroot(mg_L_D_star), mean(alpha_dist), mean(beta_dist), mean(theta_dist), ...
                               mean(R_l_dist), mean(lmb_f_dist), mean(lmb_l_dist), mg_lmb_D);
    sust_dn_mg_lmb_D = costfn_donothing(mean(beta_dist), mean(theta_dist), mean(lmb_f_dist), mean(lmb_l_dist), mean(R_l_dist));
    sust_mg_lmb_D = sust_da_mg_lmb_D ./ sust_dn_mg_lmb_D;
    
    sust_da_mg_lmb_f = costfn_doavulsion(mg_L_D_star, return_L_f_star(mg_L_D_star, mean(beta_dist), mean(R_l_dist)), L_d_star_ijk, ...
                               return_fn_squareroot(mg_L_D_star), mean(alpha_dist), mean(beta_dist), mean(theta_dist), ...
                               mean(R_l_dist), mg_lmb_f, mean(lmb_l_dist), mean(lmb_D_dist));
    sust_dn_mg_lmb_f = costfn_donothing(mean(beta_dist), mean(theta_dist), mg_lmb_f, mean(lmb_l_dist), mean(R_l_dist));
    sust_mg_lmb_f = sust_da_mg_lmb_f ./ sust_dn_mg_lmb_f;
    
    sust_da_mg_lmb_l = costfn_doavulsion(mg_L_D_star, return_L_f_star(mg_L_D_star, mean(beta_dist), mean(R_l_dist)), L_d_star_ijk, ...
                               return_fn_squareroot(mg_L_D_star), mean(alpha_dist), mean(beta_dist), mean(theta_dist), ...
                               mean(R_l_dist), mean(lmb_f_dist), mg_lmb_l, mean(lmb_D_dist));
    sust_dn_mg_lmb_l = costfn_donothing(mean(beta_dist), mean(theta_dist), mean(lmb_f_dist), mg_lmb_l, mean(R_l_dist));
    sust_mg_lmb_l = sust_da_mg_lmb_l ./ sust_dn_mg_lmb_l;
    
    sust_da_mg_alpha = costfn_doavulsion(mg_L_D_star, return_L_f_star(mg_L_D_star, mean(beta_dist), mean(R_l_dist)), L_d_star_ijk, ...
                               return_fn_squareroot(mg_L_D_star), mg_alpha, mean(beta_dist), mean(theta_dist), ...
                               mean(R_l_dist), mean(lmb_f_dist), mean(lmb_l_dist), mean(lmb_D_dist));
    sust_dn_mg_alpha = costfn_donothing(mean(beta_dist), mean(theta_dist), mean(lmb_f_dist), mean(lmb_l_dist), mean(R_l_dist));
    sust_mg_alpha = sust_da_mg_alpha ./ sust_dn_mg_alpha;
    
    sust_da_mg_beta = costfn_doavulsion(mg_L_D_star, return_L_f_star(mg_L_D_star, mg_beta, mean(R_l_dist)), L_d_star_ijk, ...
                               return_fn_squareroot(mg_L_D_star), mean(alpha_dist), mg_beta, mean(theta_dist), ...
                               mean(R_l_dist), mean(lmb_f_dist), mean(lmb_l_dist), mean(lmb_D_dist));
    sust_dn_mg_beta = costfn_donothing(mg_beta, mean(theta_dist), mean(lmb_f_dist), mean(lmb_l_dist), mean(R_l_dist));
    sust_mg_beta = sust_da_mg_beta ./ sust_dn_mg_beta;
    
    sust_da_mg_theta = costfn_doavulsion(mg_L_D_star, return_L_f_star(mg_L_D_star, mean(beta_dist), mean(R_l_dist)), L_d_star_ijk, ...
                               return_fn_squareroot(mg_L_D_star), mean(alpha_dist), mean(beta_dist), mg_theta, ...
                               mean(R_l_dist), mean(lmb_f_dist), mean(lmb_l_dist), mean(lmb_D_dist));
    sust_dn_mg_theta = costfn_donothing(mean(beta_dist), mg_theta, mean(lmb_f_dist), mean(lmb_l_dist), mean(R_l_dist));
    sust_mg_theta = sust_da_mg_theta ./ sust_dn_mg_theta;
    
    sust_da_mg_R_l = costfn_doavulsion(mg_L_D_star, return_L_f_star(mg_L_D_star, mean(beta_dist), mg_R_l), L_d_star_ijk, ...
                               return_fn_squareroot(mg_L_D_star), mean(alpha_dist), mean(beta_dist), mean(theta_dist), ...
                               mg_R_l, mean(lmb_f_dist), mean(lmb_l_dist), mean(lmb_D_dist));
    sust_dn_mg_R_l = costfn_donothing(mean(beta_dist), mean(theta_dist), mean(lmb_f_dist), mean(lmb_l_dist), mg_R_l);
    sust_mg_R_l = sust_da_mg_R_l ./ sust_dn_mg_R_l;
    
    %subplot(2,2,2); hold on;
    pspacecmap = winter(20);
    pspacecint = 0:0.05:2;
    pspaceclim = [0.3, 1.25];
    sustFIG = figure(); 
    subplot(2, 4, 1); hold on;
        plot([1, 1], [0, 1e-1], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        [C,h] = contour(mg_L_D_star, mg_lmb_D, sust_mg_lmb_D, pspacecint, 'LineWidth', 2);
        view([0 90])
        clabel(C,h, 'interpreter', label_interpreter, 'Color', [0.5 0.5 0.5])
        caxis(pspaceclim)
        text(0.85, 0.85, 'a', 'units', 'normalized', 'FontSize', 12, 'BackgroundColor', [1 1 1])
        %xlabel('${L_A}^*$', 'interpreter', label_interpreter)
        ylabel('$\lambda_D$', 'interpreter', label_interpreter)
        colormap(pspacecmap)
        axis square
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
    subplot(2, 4, 2); hold on;
        plot([1, 1], [0, 10], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        [C,h] = contour(mg_L_D_star, mg_lmb_f, sust_mg_lmb_f, pspacecint, 'LineWidth', 2);
        view([0 90])
        clabel(C,h, 'interpreter', label_interpreter, 'Color', [0.5 0.5 0.5])
        caxis(pspaceclim)
        text(0.85, 0.85, 'b', 'units', 'normalized', 'FontSize', 12, 'BackgroundColor', [1 1 1])
        %xlabel('${L_A}^*$', 'interpreter', label_interpreter)
        ylabel('$\lambda_f$', 'interpreter', label_interpreter)
        colormap(pspacecmap)
        axis square
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
    subplot(2, 4, 3); hold on;
        plot([1, 1], [0, 10], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        [C,h] = contour(mg_L_D_star, mg_lmb_l, sust_mg_lmb_l, pspacecint, 'LineWidth', 2);
        view([0 90])
        clabel(C,h, 'interpreter', label_interpreter, 'Color', [0.5 0.5 0.5])
        caxis(pspaceclim)
        text(0.85, 0.85, 'c', 'units', 'normalized', 'FontSize', 12, 'BackgroundColor', [1 1 1])
        %xlabel('${L_A}^*$', 'interpreter', label_interpreter)
        ylabel('$\lambda_l$', 'interpreter', label_interpreter)
        colormap(pspacecmap)
        axis square
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
    subplot(2, 4, 4); hold on;
        plot([1, 1], [0, 1e-1], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        [C,h] = contour(mg_L_D_star, mg_alpha, sust_mg_alpha, pspacecint, 'LineWidth', 2);
        view([0 90])
        clabel(C,h, 'interpreter', label_interpreter, 'Color', [0.5 0.5 0.5])
        caxis(pspaceclim)
        text(0.85, 0.85, 'd', 'units', 'normalized', 'FontSize', 12, 'BackgroundColor', [1 1 1])
        xlabel('${L_D}^*$', 'interpreter', label_interpreter)
        ylabel('$\alpha$', 'interpreter', label_interpreter)
        colormap(pspacecmap)
        axis square
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
    subplot(2, 4, 5); hold on;
        plot([1, 1], [0.4, 0.6], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        [C,h] = contour(mg_L_D_star, mg_beta, sust_mg_beta, pspacecint, 'LineWidth', 2);
        view([0 90])
        clabel(C,h, 'interpreter', label_interpreter, 'Color', [0.5 0.5 0.5])
        caxis(pspaceclim)
        text(0.85, 0.85, 'e', 'units', 'normalized', 'FontSize', 12, 'BackgroundColor', [1 1 1])
        xlabel('${L_D}^*$', 'interpreter', label_interpreter)
        ylabel('$\beta$', 'interpreter', label_interpreter)
        colormap(pspacecmap)
        axis square
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
    subplot(2, 4, 6); hold on;
        plot([1, 1], [0, 45], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        [C,h] = contour(mg_L_D_star, mg_theta, sust_mg_theta, pspacecint, 'LineWidth', 2);
        view([0 90])
        clabel(C,h, 'interpreter', label_interpreter, 'Color', [0.5 0.5 0.5])
        caxis(pspaceclim)
        text(0.85, 0.85, 'f', 'units', 'normalized', 'FontSize', 12, 'BackgroundColor', [1 1 1])
        xlabel('${L_D}^*$', 'interpreter', label_interpreter)
        ylabel('$\theta$', 'interpreter', label_interpreter)
        colormap(pspacecmap)
        axis square
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
    subplot(2, 4, 7); hold on;
        plot([1, 1], [0, 0.5], '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2)
        [C,h] = contour(mg_L_D_star, mg_R_l, sust_mg_R_l, pspacecint, 'LineWidth', 2);
        view([0 90])
        clabel(C,h, 'interpreter', label_interpreter, 'Color', [0.5 0.5 0.5])
        caxis(pspaceclim)
        text(0.85, 0.85, 'g', 'units', 'normalized', 'FontSize', 12, 'BackgroundColor', [1 1 1])
        xlabel('${L_D}^*$', 'interpreter', label_interpreter)
        ylabel('$R_l$', 'interpreter', label_interpreter)
        colormap(pspacecmap)
        axis square
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
    set(sustFIG, 'Pos', [700 500 1118 512]);
    set(sustFIG, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    if printOut
        print('-dpng', '-r300', [print_root, 'engineered_sustnum_contours.png']);
        set(gcf, 'InvertHardcopy', 'on')
        print('-depsc', '-opengl', '-r300', [print_root, 'engineered_sustnum_contours.eps']);
        pause(1);
    end
    
    
    if false
        figure(); hold on
            surf(mg_L_D_star, mg_lmb_f, sust_mg_2, 'EdgeColor', 'none')
            contour3(mg_L_D_star, mg_lmb_f, sust_mg_2, [1 1], 'w-', 'LineWidth', 2)
            [C,h] = contour3(mg_L_D_star, mg_lmb_f, sust_mg_2, 0.7:0.05:1.4, '-', 'LineWidth', 1, 'Color', [1 1 1]);
            clabel(C,h, 'interpreter', label_interpreter, 'Color', [1 1 1])
            view([0 90])
            %colormap(ndvi_colormap)
            caxis([0.87 1.2])
        xlabel('${L_D}^*$', 'interpreter', label_interpreter)
        ylabel('$\lambda_f$', 'interpreter', label_interpreter)
        axis square
        set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
        box on
        set(gca, 'Layer', 'top')
    end
    

    %% components of the sustnum increase with lambda_f
    sncmptFIG = figure(); hold on
        for i = 1:length(lamb_fs_plot)
            f(i) = fill([L_D_star_list, fliplr(L_D_star_list)], [mced_lmb_f_mean(:,i)+mced_lmb_f_std(:,i); flipud(mced_lmb_f_mean(:,i)-mced_lmb_f_std(:,i))], ...
                    lmbf_cmap(i,:), 'EdgeColor', factor_color(lmbf_cmap(i,:), 1.1), 'FaceAlpha', 0.2);
        end
        for i = 1:length(lamb_fs_plot)
            plot(L_D_star_list, mced_lmb_f_mean(:,i), '-', 'Color', factor_color(lmbf_cmap(i,:),0.8), 'LineWidth', 1.2)
            plot(L_D_star_list, mced_lmb_f_dn_mean(:,i), '--', 'Color', factor_color(lmbf_cmap(i,:),0.8), 'LineWidth', 1.2)
        end
    legend([f'], arrayfun(@(x) ['$\lambda_f =~$', num2str(x)], lamb_fs_plot, 'unif', 0), 'location', 'southwest', 'interpreter', 'latex')
    ylim([0.5, 1.25])
    xlim([0, 1.7])
    yticks(0.5:0.25:1.25)
    ylabel('${\lambda_\Pi}/{T_A}^*$')    
    xlabel('${L_D}^*$', 'interpreter', label_interpreter)
    axis square
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'ticklabelinterpreter', 'latex');
    box on
    set(gca, 'Layer', 'top')
    set(sncmptFIG, 'Pos', [700 500 329   256]);
    set(sncmptFIG, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'color', bgcolor)
    if printOut
        print('-dpng', '-r300', [print_root, 'engineered_sustnum_parts.png']);
        set(sncmptFIG, 'InvertHardcopy', 'on')
        print('-depsc', '-opengl', '-r300', [print_root, 'engineered_sustnum_parts.eps']);
        pause(1);
    end
    
    %% real world cost data
    yuan_to_dollars = 0.15;
    arable_km2_damage_doll = 6e6 * yuan_to_dollars;
    arable_km2_revenue_doll = 6e6 * yuan_to_dollars;

    avulsion_km_cost_doll = 2e6 * yuan_to_dollars;
    avulsion_fix_cost_doll = 180e6 * yuan_to_dollars;
    
    arable_km2_damage = arable_km2_damage_doll / arable_km2_damage_doll;
    arable_km2_revenue = arable_km2_revenue_doll / arable_km2_damage_doll;

    avulsion_km_cost = avulsion_km_cost_doll / arable_km2_damage_doll;
    avulsion_fix_cost = avulsion_fix_cost_doll / arable_km2_damage_doll;
    
    lamb_l = arable_km2_damage_doll / arable_km2_damage_doll;
    lamb_A = avulsion_fix_cost_doll / arable_km2_revenue_doll;
    
   
end

function [color_factored] = factor_color(color, factor)

    if factor > 1
        color_factored = color .* factor;
        color_factored = min(color_factored, [1 1 1]);
    elseif factor < 1
        color_factored = color .* factor;
        color_factored = max(color_factored, [0 0 0]);
    else
        color_factored = color;
    end
end