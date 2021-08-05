function goldblum_FlyDiscoCaboose_wrapper(experiment_folder_path, settings_folder_path, overriding_analysis_parameters_as_list)
    % This is the function that is submitted by goldblum to the bqueue to run each
    % experiment.
  
    % Handle arguments
    if ~exist('settings_folder_path', 'var') || isempty(settings_folder_path) ,
        script_folder_path = fileparts(mfilename('fullpath')) ;
        fly_disco_analysis_folder_path = fileparts(script_folder_path) ;
        settings_folder_path = fullfile(fly_disco_analysis_folder_path, 'settings') ;
    end
    if ~exist('overriding_analysis_parameters_as_list', 'var') || isempty(overriding_analysis_parameters_as_list) ,
        overriding_analysis_parameters_as_list = cell(1, 0) ;
    end
    
    % If get here, create the lock file
    analysis_in_progress_file_path = fullfile(experiment_folder_path, 'ANALYSIS-IN-PROGRESS') ;
    touch(analysis_in_progress_file_path) ;

    % Print the date to the stdout, so it gets logged
    dt = datetime('now') ;
    date_as_string = string(datetime(dt, 'Format', 'uuuu-MM-dd')) ;
    time_as_string = string(datetime(dt, 'Format', 'HH:mm:ss')) ;
    header_string = sprintf('Running FlyDiscoCaboose() on %s at %s', date_as_string, time_as_string) ;
    asterisks_string = repmat('*', [1 length(header_string)]) ;
    fprintf('\n\n\n\n\n') ;
    fprintf('%s\n', asterisks_string) ;    
    fprintf('%s\n', header_string) ;
    fprintf('%s\n\n', asterisks_string) ;    
    
    % Check for an already-existing failure file.  If none exists, create one.
    does_failure_file_already_exist = false ;
    analysis_errored_out_file_path = fullfile(experiment_folder_path, 'ANALYSIS-ERRORED-OUT') ;
    if exist(analysis_errored_out_file_path, 'file') ,
        does_failure_file_already_exist = true ;
    end
    analysis_incomplete_file_path = fullfile(experiment_folder_path, 'ANALYSIS-INCOMPLETE') ;
    if exist(analysis_incomplete_file_path, 'file') ,
        does_failure_file_already_exist = true ;
    end
    if does_failure_file_already_exist ,
        % do nothing
    else        
        touch(analysis_errored_out_file_path) ;
    end
    
    % Convert param list to a struct
    overriding_analysis_parameters = struct_from_name_value_list(overriding_analysis_parameters_as_list) ;

    % Build up the parameters cell array
    default_analysis_parameters = struct('settingsdir', {settings_folder_path}) ;
    
    % Combine the caller-supplied analysis parameters with the defaults       
    analysis_parameters = merge_structs(default_analysis_parameters, overriding_analysis_parameters) ;
    
    % Call the function to do the real work
    did_caboose_error_out = false ;
    try
        FlyDiscoCaboose(experiment_folder_path, analysis_parameters) ;
    catch me
        % Whatever happens, want to write out one of the two ANALYSIS-* files
        fprintf('Encountered error in FlyDiscoCaboose():\n') ;
        fprintf('%s\n', me.getReport()) ;
        did_caboose_error_out = true ;
    end
       
    % Clear the lock file
    if exist(analysis_in_progress_file_path, 'file') ,
        delete(analysis_in_progress_file_path) ;
    end
    
    % Error out if FlyDiscoCaboose() errored out or returned success==false, but
    % distinguish between them in the log.
    if did_caboose_error_out ,
        [~,experiment_folder_name] = fileparts2(experiment_folder_path) ;
        error('FlyDiscoCaboose() errored out on experiment %s!', experiment_folder_name) ;  % want to return a non-zero error code
    end        
end
