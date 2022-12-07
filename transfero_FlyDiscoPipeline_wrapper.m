function transfero_FlyDiscoPipeline_wrapper(experiment_folder_path, user_name_for_configuration_purposes, overriding_analysis_parameters_as_list, do_try)
    % This is the function that is called from transfero_FlyDiscoPipeline_wrapper_wrapper.py, 
    % which is called by Transfero once for each experiment to be analyzed.
  
    % Deal with args
    if ~exist('overriding_analysis_parameters_as_list', 'var') || isempty(overriding_analysis_parameters_as_list) ,
        overriding_analysis_parameters_as_list = cell(1, 0) ;
    end
    if ~exist('do_try', 'var') || isempty(do_try) ,
        % true means to wrap the main call to FlyDiscoPipeline in a try-catch clause
        % setting this to false is useful when debugging
        do_try = true ;
    end        

    % Load the per-lab Goldblum/FDA configuration file
    %user_name = get_user_name() ;
    configuration_function_name = sprintf('%s_configuration', user_name_for_configuration_purposes) ;
    configuration = feval(configuration_function_name) ;
    
    % Unpack the per-lab configuration file
    cluster_billing_account_name = configuration.cluster_billing_account_name ;
    host_name_from_rig_index = configuration.host_name_from_rig_index ;
    rig_user_name_from_rig_index = configuration.rig_user_name_from_rig_index ;
    data_folder_path_from_rig_index = configuration.data_folder_path_from_rig_index ;
    destination_folder = configuration.destination_folder ;    
    settings_folder_path = configuration.settings_folder_path
    %does_use_per_user_folders = configuration.does_use_per_user_folders ;
    to_process_folder_name = 'to-process' ;

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
    
    % Convert param list to a struct
    overriding_analysis_parameters = struct_from_name_value_list(overriding_analysis_parameters_as_list) ;

    % Build up the parameters cell array
    default_analysis_parameters = struct('settingsdir', {settings_folder_path}) ;
    
    % Combine the caller-supplied analysis parameters with the defaults       
    analysis_parameters_with_overrides = merge_structs(default_analysis_parameters, overriding_analysis_parameters) ;

    % Now turn off the auto-checks-complete, we do that separately, afterwards
    analysis_parameters = analysis_parameters_with_overrides ;
    analysis_parameters.doautomaticcheckscomplete = 'off' ;
    
    % Call the function to do the real work
    main_pipeline_exception_maybe = [] ;
    if do_try ,
        try
            FlyDiscoPipeline(experiment_folder_path, analysis_parameters) ;
        catch pipeline_exception ,
            main_pipeline_exception_maybe = pipeline_exception ;
        end
    else
        FlyDiscoPipeline(experiment_folder_path, analysis_parameters) ;
    end

    % Now we do the "caboose" phase, which we want to run even if the 
    % main phase fails 

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
    
    % Now turn off everything *except* the auto-checks-complete
    caboose_analysis_parameters = analysis_parameters_with_overrides ;
    caboose_analysis_parameters.doautomaticchecksincoming = 'off' ;
    caboose_analysis_parameters.doflytracking = 'off' ;
    caboose_analysis_parameters.doregistration = 'off' ;
    caboose_analysis_parameters.doledonoffdetection = 'off' ;
    caboose_analysis_parameters.dosexclassification = 'off' ;
    caboose_analysis_parameters.dotrackwings = 'off' ;
    caboose_analysis_parameters.docomputeperframefeatures = 'off' ;
    caboose_analysis_parameters.docomputehoghofperframefeatures = 'off' ;
    caboose_analysis_parameters.dojaabadetect = 'off' ;
    caboose_analysis_parameters.docomputeperframestats = 'off' ;
    caboose_analysis_parameters.doplotperframestats = 'off' ;
    caboose_analysis_parameters.domakectraxresultsmovie = 'off' ;
    caboose_analysis_parameters.doextradiagnostics = 'off' ;
    caboose_analysis_parameters.doanalysisprotocol = 'off' ;
    
    % Call the function to do the real work
    FlyDiscoPipeline(experiment_folder_path, caboose_analysis_parameters) ;
    
    % If there was an exception thrown during the main pipeline, rethrow it so that
    % we exit with an error return code when running in batch mode.
    if ~isempty(main_pipeline_exception_maybe) ,
        main_pipeline_exception = main_pipeline_exception_maybe(1) ;
        rethrow(main_pipeline_exception) ;
    end
end
