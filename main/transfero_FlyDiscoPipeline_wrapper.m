function transfero_FlyDiscoPipeline_wrapper(experiment_folder_path, user_name_for_configuration_purposes, varargin)
% This is the function that is called from transfero_FlyDiscoPipeline_wrapper_wrapper.py,
% which is called by Transfero once for each experiment to be analyzed.

% Collect the optional argument name, value pairs into a struct, for easier
% handling.
optional_arguments_as_list = varargin ;
argument_analysis_parameters = struct_from_name_value_list(optional_arguments_as_list) ;

% Load the per-lab FDA configuration file
%user_name = get_user_name() ;
configuration_function_name = sprintf('%s_configuration', user_name_for_configuration_purposes) ;
configuration = feval(configuration_function_name) ;

% Unpack the per-lab configuration file
cluster_billing_account_name = configuration.cluster_billing_account_name ;
%host_name_from_rig_index = configuration.host_name_from_rig_index ;
%rig_user_name_from_rig_index = configuration.rig_user_name_from_rig_index ;
%data_folder_path_from_rig_index = configuration.data_folder_path_from_rig_index ;
%destination_folder = configuration.destination_folder ;
settings_folder_path = configuration.settings_folder_path ;
%does_use_per_user_folders = configuration.does_use_per_user_folders ;
%to_process_folder_name = 'to-process' ;

% Print the date to the stdout, so it gets logged
dt = datetime('now') ;
date_as_string = string(datetime(dt, 'Format', 'uuuu-MM-dd')) ;
time_as_string = string(datetime(dt, 'Format', 'HH:mm:ss')) ;
header_string = sprintf('Running FlyDiscoPipeline() on %s at %s', date_as_string, time_as_string) ;
asterisks_string = repmat('*', [1 length(header_string)]) ;
fprintf('\n\n\n\n\n') ;
fprintf('%s\n', asterisks_string) ;
fprintf('%s\n', header_string) ;
fprintf('%s\n\n', asterisks_string) ;

% Build up the parameters cell array
default_analysis_parameters = struct('settingsdir', {settings_folder_path}, ...
                                     'cluster_billing_account_name', cluster_billing_account_name) ;

% Combine the caller-supplied analysis parameters with the defaults
main_phase_analysis_parameters = merge_structs(default_analysis_parameters, argument_analysis_parameters) ;

% Convert the analysis arguments back to a name-value list
main_phase_analysis_parameters_as_list = name_value_list_from_struct(main_phase_analysis_parameters) ;

% Call the function to do the real work
FlyDiscoPipeline(experiment_folder_path, main_phase_analysis_parameters_as_list{:}) ;
