function fly_disco_analysis_pipeline_wrapper(experiment_folder_path, settings_folder_path, overriding_analysis_parameters_as_list)
    % Handle arguments
    if ~exist('settings_folder_path', 'var') || isempty(settings_folder_path) ,
        script_folder_path = fileparts(mfilename('fullpath')) ;
        fly_disco_analysis_folder_path = fileparts(script_folder_path) ;
        settings_folder_path = fullfile(fly_disco_analysis_folder_path, 'settings') ;
    end
    if ~exist('overriding_analysis_parameters', 'var') || isempty(overriding_analysis_parameters_as_list) ,
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
    
%     % Read the experiment metadata to determine the analysis_protoocol
%     metadata_file_path = determine_metadata_file_path(experiment_folder_path) ;
%     analysis_protocol_folder_name = analysis_protocol_from_metadata_file(metadata_file_path, settings_folder_path) ;    
%     fprintf('Analysis protocol is: %s\n', analysis_protocol_folder_name) ;
%     analysis_protocol_folder_path = fullfile(settings_folder_path, analysis_protocol_folder_name) ;
%     escaped_analysis_protocol_folder_path = escape_path_for_bash(analysis_protocol_folder_path) ;
%     command_line = sprintf('realpath %s', escaped_analysis_protocol_folder_path) ;
%     stdout = system_with_error_handling(command_line) ;
%     canonical_analysis_protocol_folder_path = strtrim(stdout) ;
%     fprintf('Canonical path to analysis protocol folder is:\n  %s\n', canonical_analysis_protocol_folder_path) ;

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
