function fly_disco_analysis_pipeline_wrapper(experiment_folder_path, settings_folder_path, overriding_analysis_parameters_as_list)
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

    % Check for the lock file
    analysis_in_progress_file_path = fullfile(experiment_folder_path, 'ANALYSIS-IN-PROGRESS') ;
    if exist(analysis_in_progress_file_path, 'file') ,
        error('Not going to run FlyDiscoPipeline() on experiment folder:\n\n  %s\n\nANALYSIS-IN-PROGRESS file is already present.\n', ...
              experiment_folder_path) ;  % error() to return a non-zero error code
    end
    
    % If get here, create the lock file
    touch(analysis_in_progress_file_path) ;

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
    
    % Delete any pre-existing success/failure files
    analysis_successful_file_path = fullfile(experiment_folder_path, 'ANALYSIS-COMPLETED-SUCCESSFULLY') ;
    if exist(analysis_successful_file_path, 'file') ,
        delete(analysis_successful_file_path) ;
    end
    analysis_failed_file_path = fullfile(experiment_folder_path, 'ANALYSIS-FAILED') ;
    if exist(analysis_failed_file_path, 'file') ,
        delete(analysis_failed_file_path) ;
    end
    
    % Convert param list to a struct
    overriding_analysis_parameters = struct_from_name_value_list(overriding_analysis_parameters_as_list) ;

    % Build up the parameters cell array
    default_analysis_parameters = struct('settingsdir', {settings_folder_path}) ;
    
    % Combine the caller-supplied analysis parameters with the defaults       
    analysis_parameters = merge_structs(default_analysis_parameters, overriding_analysis_parameters) ;
    
    % Call the function to do the real work
    try
        [success, msgs, stage] = FlyDiscoPipeline(experiment_folder_path, analysis_parameters) ;
    catch me
        % Whatever happens, want to write out one of the two ANALYSIS-* files
        fprintf('Encountered error in FlyDiscoPipeline():\n') ;
        fprintf('%s\n', me.getReport()) ;
        success = false ;
        msgs = {} ;
        stage = '<unknown>' ;
    end
    
    % Deal with success or failure
    if success ,
        analysis_successful_file_path = fullfile(experiment_folder_path, 'ANALYSIS-COMPLETED-SUCCESSFULLY') ;
        touch(analysis_successful_file_path) ;
    else        
        analysis_failed_file_path = fullfile(experiment_folder_path, 'ANALYSIS-FAILED') ;
        touch(analysis_failed_file_path) ;
        fprintf('FlyDiscoPipeline() encountered one or more problems at stage %s:\n', stage)
        for i = 1 : length(msgs) ,
            msg = msgs{i} ;
            fprintf('%s\n', msg) ;
        end
    end
    
    % Clear the lock file
    if exist(analysis_in_progress_file_path, 'file') ,
        delete(analysis_in_progress_file_path) ;
    end
    
    % Error out if failed
    if ~success ,
        error('FlyDiscoPipeline() failed!') ;  % want to return a non-zero error code
    end        
end