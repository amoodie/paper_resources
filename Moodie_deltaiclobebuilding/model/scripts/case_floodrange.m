%% script for running the flood range scenario

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
cfg.run_n = 24; % should be 24 for real runs

% setup the avulsion parameters
cfg.Int = 1; % intermittancy factor, 1 for full time, variable discharge simulation
cfg.preavul_switch = 'any-setup+over+topset'; % avulsion setup/trigger method
cfg.postavul_switch = 'taper-to-bankfull'; % avulsion effect on the bed
cfg.preavul_thresh = 0.5;
cfg.preavul_trigg = NaN; % how much over-floodplain flow is needed before avulsion

% setup model input / boundary conditions
cfg.initial_channel = 'new'; % use new channel params or load an old set 
cfg.mouth_switch = 'bartip'; % mouth progradation method
cfg.Qw_str = 'engineered'; % how to determine the discharge scheme

% setup saving/output params
cfg.storage = false;
cfg.savedata_filename = NaN; % fullfile('..', 'output_data', 'longterm_data.mat');
cfg.savefinal_filename = NaN;
cfg.analysis_filename = fullfile('..', 'output_data', 'floodrange_analysis_lastonly.mat');

% setup runtime information outputs
cfg.runtime_plot = false;
cfg.runtime_save = false;
cfg.runtime_plot_int = 30;
cfg.runtime_stdout_int = 100;

flood_list = [6000]; % [500, 1000:1000:6000]; % the peak discharge modeled during the flood

% run the model with cfg
floodrange = cell(length(flood_list), 1);
for i = 1:length(flood_list)
    
    cfg.Qw_num = flood_list(i); % numeric argument to Qw func
    
    [~, analysis] = virtualdelta(cfg);
    floodrange(i) = {analysis};
    disp(['floodrange case = ' num2str(i) ' of ' num2str(length(flood_list))])

end

save('../output_data/floodrange_analysis.mat', 'floodrange', 'flood_list')
