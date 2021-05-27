%% script for running the long-term case scenario

% add things to path
source_path = genpath(fullfile('..', 'source'));
input_data_path = genpath(fullfile('..', 'input_data'));
post_proc_path = genpath(fullfile('..', 'post_processing'));
addpath(source_path, input_data_path, post_proc_path)


%% model configuration params
% initialize Yellow River configuration parameters
[cfg] = yellowriver_configuration();

% setup how to run the timing scheme
cfg.run_str = 'nAvuls';

% setup the avulsion parameters
cfg.preavul_switch = 'any-setup+topset'; % avulsion setup/trigger method
cfg.postavul_switch = 'taper-to-bankfull'; % avulsion effect on the bed
cfg.preavul_trigg = NaN; % how much over-floodplain flow is needed before avulsion

% setup model input / boundary conditions
cfg.Int = 1; % intermittency factor
cfg.initial_channel = 'new'; % use new channel params or load an old set 
cfg.mouth_switch = 'bartip'; % mouth progradation method
cfg.Qw_str = 'mean_historical'; % how to determine the discharge scheme
cfg.Qw_num = NaN; % numeric argument to Qw func

% setup saving/output params
cfg.storage = true;

% setup runtime information outputs
cfg.runtime_plot = false;
cfg.runtime_save = false;
cfg.runtime_plot_int = 30;
cfg.runtime_stdout_int = 100;

%% initial setup runs here
spin_catalog = [0.4, 14; 0.4, 15; 0.45, 14; 0.45, 15; 0.5, 14; 0.5, 15];

if false
    for i = 1:size(spin_catalog, 1)
        cfg.preavul_thresh = spin_catalog(i,1);
        cfg.run_n = spin_catalog(i,2);
        
        cfg.savedata_filename = fullfile('..', 'output_data', 'engineered', 'spinups', ['engineered_spinup_data_single_hist_', num2str(i), '.mat']);
        cfg.savefinal_filename = fullfile('..', 'output_data', 'engineered', 'spinups', ['engineered_spinup_data_single_', num2str(i), '.mat']);
        cfg.analysis_filename = fullfile('..', 'output_data', 'engineered', 'spinups', ['engineered_spinup_analysis_single_', num2str(i), '.mat']);

        [~, ~] = virtualdelta(cfg);
    end
end

%% setup and run the second part off the base of the above run
cfg.run_n = 1; % only one avulsion follows the enginneered avulsion
cfg.initial_channel = 'engineered'; % use new channel params or load an old set
cfg.savedata_filename = NaN;
cfg.savefinal_filename = NaN;
cfg.analysis_filename = fullfile('..', 'output_data', 'engineered', 'engineered_analysis_secondpart.mat');

avulsion_location_list = [0.05:0.05:0.5, 0.6:0.1:1.6]; % location of where to trigger the eng avulsion
n_expts = 10; % number of trials per avulsion location
save(fullfile('..', 'output_data', 'engineered', 'engineered_single_avulsion_location_list.mat'), 'avulsion_location_list', 'n_expts')

% run the model with cfg
avulsion_location_range = cell(length(avulsion_location_list), n_expts);
for i = 1:length(avulsion_location_list)

    % grab new avulsion location (LA)
    avulsion_location = avulsion_location_list(i);
    
    for j = 1:n_expts
        % repeat the experiment n_expts times, varying the avulsion chars by
        %   a small amount to get some uncertainty sense
        
        % load the old spinup output
        spin_rand = randi(size(spin_catalog,1));
        eng = load(fullfile('..', 'output_data', 'engineered', 'spinups', ['engineered_spinup_data_single_', num2str(spin_rand), '.mat']));
        
        % calculate avulloc for new LA
        LA = eng.Lblow * avulsion_location; % length of avulsion, m
        avul_loc_m = (eng.mou0_idx*eng.dx) - LA;
        rand_offset = 0; % % randi(3)-2;
        eng.avulloc = get_idx(eng.x, avul_loc_m) + rand_offset; % index of avulsion length to use
        diversion_index = eng.avulloc;

        % do avulsion on the delta with LA
            % complete avulsion in planform domain
            % ALREADY HAPPENED IN SAVED FILE!!

            % complete avulsion in long profile
            % eng.eta0 = eng.s.eta(:,end-1); % don't need to do, eta0 is already profile just before avulsion
            random_change = 2*( randi(3)-2 );
            [eng.eta] = avulsion_long(eng.eta_i, eng.eta0, eng.zed, eng.rad_idx, eng.Hnbf, eng.avulloc, eng.dx, eng.cfg.postavul_switch, random_change);
            if avulsion_location < 0.5
                lobe_bed = (eng.zed - eng.Hnbf);
                % eng.eta(eng.avulloc:eng.mou0_idx) = lobe_bed(eng.avulloc:eng.mou0_idx);
                eng.eta(eng.avulloc:(eng.avulloc+floor(eng.Bo0/eng.dx))) = lobe_bed(eng.avulloc:(eng.avulloc+floor(eng.Bo0/eng.dx)));
            end
            [eng.Bc] = set_Bc(eng.mou_idx, eng.Bc0, eng.thet, eng.nx, eng.dx);
            [eng.Be] = set_Be(eng.Bc0, eng.Bf0, eng.Bo0, eng.mou_idx, eng.rad_idx, eng.nx);
            idx_upstream_avulsion = find(eng.eta0 ~= eng.eta, 1, 'first');

            % reset/update params
            eng.lastavul_idx = eng.cnt.delap;
            eng.Hbasin = eng.H0 - eng.eta;
            eng.cnt.nAvuls = eng.cnt.nAvuls + 1;
            eng.eta_lastavul = eng.eta;
            eng.flux_to_lobe = 0;
            eng.avulsion_this_year = true;

            % save it for loading in the next step in virtualdelta
            eta = eng.eta;
            zed = eng.zed;
            rad = eng.rad;
            save(fullfile('..', 'output_data', 'engineered', 'engineered_input_lastonly.mat'), ...
                'eta', 'zed', 'rad')

            % clear some memory
            clear eng eta zed rad
            pause(2)

        cfg.preavul_thresh = spin_catalog(spin_rand,1);
        [data, analysis] = virtualdelta(cfg);
        data.idx_upstream_avulsion = idx_upstream_avulsion;
        data.diversion_idx = diversion_index;
        %avulsion_location_range(i,j) = {{data, analysis}};

        save(fullfile('..', 'output_data', 'engineered', ['engineered_data_single_LA',num2str(i),'_expt',num2str(j),'.mat']), ...
            'data')
    end

    disp(['avulsion_location_range case = ' num2str(i) ' of ' num2str(length(avulsion_location_list))])

    clear data analysis
end

save(fullfile('..', 'output_data', 'engineered', 'engineered_single_data_analysis.mat'), ...
    'avulsion_location_range', 'avulsion_location_list', '-v7.3')
