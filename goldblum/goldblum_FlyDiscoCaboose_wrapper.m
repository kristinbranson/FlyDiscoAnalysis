function goldblum_FlyDiscoCaboose_wrapper(experiment_folder_path, settings_folder_path, overriding_analysis_parameters_as_list)
    % This is the function that is submitted by goldblum to the bqueue to run each
    % experiment.
  
    % Handle arguments
    if ~exist('settings_folder_path', 'var') || isempty(settings_folder_path) ,
        settings_folder_path = default_settings_folder_path() ;
    end
    if ~exist('overriding_analysis_parameters_as_list', 'var') || isempty(overriding_analysis_parameters_as_list) ,
        overriding_analysis_parameters_as_list = cell(1, 0) ;
    end
    
    % Print the date to the stdout, so it gets logged
    dt = datetime('now') ;
    date_as_string = string(datetime(dt, 'Format', 'uuuu-MM-dd')) ;
    time_as_string = string(datetime(dt, 'Format', 'HH:mm:ss')) ;
    header_string = sprintf('Running FlyDiscoPipeline() with caboose-appropriate settings on %s at %s', date_as_string, time_as_string) ;
    asterisks_string = repmat('*', [1 length(header_string)]) ;
    fprintf('\n\n\n\n\n') ;
    fprintf('%s\n', asterisks_string) ;    
    fprintf('%s\n', header_string) ;
    fprintf('%s\n\n', asterisks_string) ;    
    
    % Convert param list to a struct
    overriding_analysis_parameters = struct_from_name_value_list(overriding_analysis_parameters_as_list) ;

    % Build up the parameters cell array
    default_analysis_parameters = struct('settingsdir', {settings_folder_path}) ;
    
    % Combine the caller-supplied analysis parameters with the defaults       
    analysis_parameters = merge_structs(default_analysis_parameters, overriding_analysis_parameters) ;
    
    % Now turn off everything *except* the auto-checks-complete
    analysis_parameters.doautomaticchecksincoming = 'off' ;
    analysis_parameters.doflytracking = 'off' ;
    analysis_parameters.doregistration = 'off' ;
    analysis_parameters.doledonoffdetection = 'off' ;
    analysis_parameters.dosexclassification = 'off' ;
    analysis_parameters.dotrackwings = 'off' ;
    analysis_parameters.docomputeperframefeatures = 'off' ;
    analysis_parameters.docomputehoghofperframefeatures = 'off' ;
    analysis_parameters.dojaabadetect = 'off' ;
    analysis_parameters.docomputeperframestats = 'off' ;
    analysis_parameters.doplotperframestats = 'off' ;
    analysis_parameters.domakectraxresultsmovie = 'off' ;
    analysis_parameters.doextradiagnostics = 'off' ;
    analysis_parameters.doanalysisprotocol = 'off' ;
    
    % Call the function to do the real work
    FlyDiscoPipeline(experiment_folder_path, analysis_parameters) ;
end
