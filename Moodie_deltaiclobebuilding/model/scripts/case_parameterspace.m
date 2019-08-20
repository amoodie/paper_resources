%% script for running the map of parameter space scenario

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
cfg.preavul_switch = 'any-setup+topset'; % avulsion setup/trigger method
cfg.postavul_switch = 'taper-to-bankfull'; % avulsion effect on the bed
cfg.preavul_trigg = NaN;

% setup model input / boundary conditions
cfg.initial_channel = 'new'; % use new channel params or load an old set 
cfg.mouth_switch = 'bartip'; % mouth progradation method

% setup saving/output params
cfg.storage = false;
cfg.savedata_filename = NaN; % fullfile('..', 'output_data', 'longterm_data.mat');
cfg.savefinal_filename = NaN;
cfg.analysis_filename = fullfile('..', 'output_data', 'paramspace_analysis_lastonly.mat');

% setup runtime information outputs
cfg.runtime_plot = false;
cfg.runtime_save = false;
cfg.runtime_plot_int = 30;
cfg.runtime_stdout_int = 100;

% set up value matrix
preavul_thresh_list = [0.2, 0.3, 0.4, 0.5, 0.6]; % how much aggradation is required before avulsion
discharge_list = [NaN, 2000, 2500, 3000, 3500]; % discharge to use, NaN is for mean run, nums is for setting flood
[paramspace_setup, paramspace_discharge] = meshgrid(preavul_thresh_list, discharge_list);
paramspace = cell(length(discharge_list), length(preavul_thresh_list));

% run row 1 of the paramspace, where discharge is mean_historical
cfg.Qw_str = 'mean_historical'; % how to determine the discharge scheme
cfg.Qw_num = NaN; % numeric argument to Qw func
for i = 1:length(preavul_thresh_list)
    
    disp(['started paramspace case = ' num2str(i) ' of ' num2str(numel(paramspace_setup))])
    cfg.preavul_thresh = paramspace_setup(1, i);
    
    [~, analysis] = virtualdelta(cfg);
%     analysis  = {paramspace_setup(1, i), NaN};
    paramspace(1, i) = {analysis};
    disp(['finished paramspace case = ' num2str(i) ' of ' num2str(numel(paramspace_setup))])
    
end

% run remaining rows of paramspace where discharge is variable
cfg.Qw_str = 'engineered'; % how to determine the discharge scheme
cfg.Qw_num = NaN; % numeric argument to Qw func
cellidxs = sub2ind(size(paramspace_setup), repelem([2:size(paramspace_setup, 1)], size(paramspace_setup, 2))', ...
                                           repmat([1:size(paramspace_setup, 2)]', size(paramspace_setup, 1)-1, 1));
                                       
cellidxs = 1:length(preavul_thresh_list);
                                       
for i = 1:length(cellidxs)
    
    disp(['started paramspace case = ' num2str(i+size(paramspace_setup,2)) ' of ' num2str(numel(paramspace_setup))])
    
    cellidx = cellidxs(i);
    
    cfg.preavul_thresh = paramspace_setup(cellidx);
    cfg.Qw_num = paramspace_discharge(cellidx);
    
    [~, analysis] = virtualdelta(cfg);
    paramspace(cellidx) = {analysis};
    disp(['finished paramspace case = ' num2str(i+size(paramspace_setup,2)) ' of ' num2str(numel(paramspace_setup))])
    
end

save('../output_data/paramspace_analysis_2500run.mat', 'paramspace', ...
                                               'preavul_thresh_list', 'discharge_list', ...
                                               'paramspace_setup', 'paramspace_discharge')
