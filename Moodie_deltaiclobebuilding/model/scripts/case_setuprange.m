%% script for running the setup range scenario

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
cfg.preavul_switch = 'any-setup+topset'; % avulsion setup/trigger method
cfg.postavul_switch = 'taper-to-bankfull'; % avulsion effect on the bed
cfg.preavul_trigg = NaN; % how much over-floodplain flow is needed before avulsion

% setup model input / boundary conditions
cfg.Int = 1; % intermittancy factor, 1 for full time, variable discharge simulation
cfg.initial_channel = 'new'; % use new channel params or load an old set
cfg.mouth_switch = 'bartip'; % mouth progradation method
cfg.Qw_str = 'mean_historical'; % how to determine the discharge scheme
cfg.Qw_num = NaN; % numeric argument to Qw func

% setup saving/output params
cfg.storage = false;
cfg.savedata_filename = NaN; % fullfile('..', 'output_data', 'setuprange_data.mat');
cfg.savefinal_filename = NaN;
cfg.analysis_filename = fullfile('..', 'output_data', 'setuprange_analysis_lastonly.mat');

% setup runtime information outputs
cfg.runtime_plot = false;
cfg.runtime_save = false;
cfg.runtime_plot_int = 30;
cfg.runtime_stdout_int = 100;

preavul_thresh_list = [0.1:0.1:1.0]; % how much aggradation is required before avulsion

% run the model with cfg
setuprange = cell(length(preavul_thresh_list), 1);
for i = 1:length(preavul_thresh_list)
    
    cfg.preavul_thresh = preavul_thresh_list(i);
    
    [~, analysis] = virtualdelta(cfg);
    setuprange(i) = {analysis};
    disp(['setuprange case = ' num2str(i) ' of ' num2str(length(preavul_thresh_list))])
    
end

save('../output_data/setuprange_analysis.mat', 'setuprange', 'preavul_thresh_list')
