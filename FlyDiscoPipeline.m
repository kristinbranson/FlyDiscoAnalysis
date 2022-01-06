function FlyDiscoPipeline(expdir, varargin)
    % Runs the FlyDisco pipeline.
    
    % Get info about the state of the repo, output to stdout
    this_script_path = mfilename('fullpath') ;
    source_folder_path = fileparts(this_script_path) ;
    git_report = get_git_report(source_folder_path) ;
    fprintf('%s', git_report) ;
    
    % Process varargin parameters
    if length(varargin)==1 && isstruct(varargin{1}) ,
        argument_parameters = varargin{1} ;
    else
        argument_parameters = struct_from_name_value_list(varargin) ;
    end
    
    % First, get the settingsdir
    settingsdir = lookup_in_struct(argument_parameters, 'settingsdir', internal_settings_folder_path()) ;
    
    % Get the metadata file path
    metadata_file_path = determine_metadata_file_path(expdir) ;
    
    % See if analysis protocol was passed in
    try
        analysis_protocol = argument_parameters.analysis_protocol ;
    catch me
        if isequal(me.identifier, 'MATLAB:nonExistentField') ,
            % Read the experiment metadata to determine the analysis_protoocol
            % Also depends on what options exist in the settingsdir
            % Don't want to do this if analysis-protocol was specified in varargin,
            % because sometimes metadata doesn't specify analysis_protocol at all.
            analysis_protocol = analysis_protocol_from_metadata_file(metadata_file_path, settingsdir) ;
        else
            rethrow(me) ;
        end
    end
    
    % Read in the analysis protocol parameters from the analysis-protocol folder, if
    % it exists
    analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
    analysis_protocol_parameters_file_path = fullfile(analysis_protocol_folder_path, 'analysis-protocol-parameters.txt') ;
    if exist(analysis_protocol_parameters_file_path, 'file') ,
        analysis_protocol_parameters = ReadParams(analysis_protocol_parameters_file_path) ;  % a struct
    else
        analysis_protocol_parameters = struct() ;
    end
    
    % Set the default analysis parameters
    default_analysis_parameters_as_list = FlyDiscoPipelineDefaultAnalysisParameters() ;
    
    %put back in when hist is ready
    %      'requiredfiles_computeperframestats',{'statsperframetxtfilestr','statsperframematfilestr','histperframetxtfilestr','histperframematfilestr'},...
    default_analysis_parameters = struct_from_name_value_list(default_analysis_parameters_as_list) ;
    
    % Combine the default parameters with those from the analysis-protocol folder and those in the arguments
    % Precedence is: argument_parameters > analysis-protocol paramters > default parameters
    analysis_parameters = merge_structs(default_analysis_parameters, analysis_protocol_parameters, argument_parameters) ;
    
    % Assign the paramters to individual variables
    datalocparamsfilestr = lookup_in_struct(analysis_parameters, 'datalocparamsfilestr') ;
    automaticchecksincoming_params = lookup_in_struct(analysis_parameters, 'automaticchecksincoming_params') ;
    registration_params = lookup_in_struct(analysis_parameters, 'registration_params') ;
    sexclassification_params = lookup_in_struct(analysis_parameters, 'sexclassification_params') ;
    computeperframefeatures_params = lookup_in_struct(analysis_parameters, 'computeperframefeatures_params') ;
    plotperframestats_params = lookup_in_struct(analysis_parameters, 'plotperframestats_params') ;    
    makectraxresultsmovie_params = lookup_in_struct(analysis_parameters, 'makectraxresultsmovie_params') ;
    automaticcheckscomplete_params = lookup_in_struct(analysis_parameters, 'automaticcheckscomplete_params') ;
    doautomaticchecksincoming = lookup_in_struct(analysis_parameters, 'doautomaticchecksincoming') ;
    doflytracking = lookup_in_struct(analysis_parameters, 'doflytracking') ;
    doregistration = lookup_in_struct(analysis_parameters, 'doregistration') ;
    doledonoffdetection = lookup_in_struct(analysis_parameters, 'doledonoffdetection') ;
    dosexclassification = lookup_in_struct(analysis_parameters, 'dosexclassification') ;
    %dotrackwings = lookup_in_struct(analysis_parameters, 'dotrackwings') ;
    docomputeperframefeatures = lookup_in_struct(analysis_parameters, 'docomputeperframefeatures') ;
    docomputehoghofperframefeatures = lookup_in_struct(analysis_parameters, 'docomputehoghofperframefeatures') ;
    dojaabadetect = lookup_in_struct(analysis_parameters, 'dojaabadetect') ;
    docomputeperframestats = lookup_in_struct(analysis_parameters, 'docomputeperframestats') ;
    doplotperframestats = lookup_in_struct(analysis_parameters, 'doplotperframestats') ;
    domakectraxresultsmovie = lookup_in_struct(analysis_parameters, 'domakectraxresultsmovie') ;
    % doextradiagnostics = lookup_in_struct(analysis_parameters, 'doextradiagnostics') ;
    % doanalysisprotocol = lookup_in_struct(analysis_parameters, 'doanalysisprotocol') ;
    doautomaticcheckscomplete = lookup_in_struct(analysis_parameters, 'doautomaticcheckscomplete') ;
    requiredfiles_automaticchecksincoming = lookup_in_struct(analysis_parameters, 'requiredfiles_automaticchecksincoming') ;
    requiredfiles_flytracker = lookup_in_struct(analysis_parameters, 'requiredfiles_flytracker') ;
    requiredfiles_registration = lookup_in_struct(analysis_parameters, 'requiredfiles_registration') ;
    requiredfiles_ledonoffdetection = lookup_in_struct(analysis_parameters, 'requiredfiles_ledonoffdetection') ;
    requiredfiles_sexclassification = lookup_in_struct(analysis_parameters, 'requiredfiles_sexclassification') ;
    %requiredfiles_wingtracking = lookup_in_struct(analysis_parameters, 'requiredfiles_wingtracking') ;
    requiredfiles_computeperframefeatures = lookup_in_struct(analysis_parameters, 'requiredfiles_computeperframefeatures') ;
    requiredfiles_computeperframestats = lookup_in_struct(analysis_parameters, 'requiredfiles_computeperframestats') ;
    requiredfiles_computehoghofperframefeatures = lookup_in_struct(analysis_parameters, 'requiredfiles_computehoghofperframefeatures') ;
    requiredfiles_plotperframestats = lookup_in_struct(analysis_parameters, 'requiredfiles_plotperframestats') ;    
    requiredfiles_makectraxresultsmovie = lookup_in_struct(analysis_parameters, 'requiredfiles_makectraxresultsmovie') ;
    requiredfiles_automaticcheckscomplete = lookup_in_struct(analysis_parameters, 'requiredfiles_automaticcheckscomplete') ;
    
    % Read in the dataloc params
    datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
    dataloc_params = ReadParams(datalocparamsfile);
    
    % Coerse all the do* variables to on/off/force
    doautomaticchecksincoming = coerce_to_on_off_force(doautomaticchecksincoming) ;
    doflytracking = coerce_to_on_off_force(doflytracking) ;
    doregistration = coerce_to_on_off_force(doregistration) ;
    doledonoffdetection = coerce_to_on_off_force(doledonoffdetection) ;
    dosexclassification = coerce_to_on_off_force(dosexclassification) ;
    %dotrackwings = coerce_to_on_off_force(dotrackwings) ;
    docomputeperframefeatures = coerce_to_on_off_force(docomputeperframefeatures) ;
    docomputehoghofperframefeatures = coerce_to_on_off_force(docomputehoghofperframefeatures) ;
    dojaabadetect = coerce_to_on_off_force(dojaabadetect) ;
    docomputeperframestats = coerce_to_on_off_force(docomputeperframestats) ;
    doplotperframestats = coerce_to_on_off_force(doplotperframestats) ; 
    domakectraxresultsmovie = coerce_to_on_off_force(domakectraxresultsmovie) ;
    % doextradiagnostics = coerce_to_on_off_force(doextradiagnostics) ; %#ok<NASGU>
    % doanalysisprotocol = coerce_to_on_off_force(doanalysisprotocol) ; %#ok<NASGU>
    doautomaticcheckscomplete = coerce_to_on_off_force(doautomaticcheckscomplete) ;
    
    % Print the settings values in use
    variable_names_to_print = ...
        { 'settingsdir', ...
        'analysis_protocol', ...
        'doautomaticchecksincoming', ...
        'doflytracking',  ...
        'doregistration', ...
        'doledonoffdetection', ...
        'dosexclassification', ...
        'docomputeperframefeatures', ...
        'docomputehoghofperframefeatures', ...
        'dojaabadetect', ...
        'docomputeperframestats', ...
        'doplotperframestats', ...
        'domakectraxresultsmovie', ...
        'doautomaticcheckscomplete' }' ;
%        'doextradiagnostics', ...
%        'doanalysisprotocol', ...
    fprintf('Settings values in FlyDiscoPipeline():\n') ;
    for i = 1 : length(variable_names_to_print) ,
        variable_name = variable_names_to_print{i} ;
        value = eval(variable_name) ;
        fprintf('  %s: %s\n', variable_name, char_array_from_value(value)) ;
    end
    fprintf('\n') ;
    
    % Print the canonical path to the analysis folder
    canonical_analysis_protocol_folder_path = realpath(absolute_path(analysis_protocol_folder_path)) ;
    fprintf('Canonical path to analysis protocol folder is:\n  %s\n\n', canonical_analysis_protocol_folder_path) ;
    
    % Print the canonical path to the experiment folder
    canonical_experiment_folder_path = realpath(absolute_path(expdir)) ;
    fprintf('Canonical path to experiment folder is:\n  %s\n\n', canonical_experiment_folder_path) ;
    
    %% check that experiment exists
    if ~exist(expdir, 'dir') ,
        error('Experiment directory %s does not exist', expdir) ;
    end
    
    
    
    % Fundamentally, the job of FlyDiscoPipeline() is to generate all the
    % appropriate analysis files from the `raw' data files. If it generates all the
    % appropriate analysis files, that's a successful run. So if the auto-checks
    % incoming (ACI) code finds that there are dead flies, then the appropriate
    % files to generate for that experiment are just
    % automatic_checks_incoming_results.txt, automatic_checks_incoming_info.mat,
    % automatic_checks_complete_results.txt, and automatic_checks_complete_info.mat.
    % If that is done successfully, then FlyDiscoPipline() should return normally,
    % and not throw an error. Thus the throwing of uncaught errors should be
    % restricted to cases where it's not even clear to the programmer how to proceed
    % to generate the output files for a particular stage. The simplest case would
    % be if a write to disk of one of the output files fails, that should clearly
    % lead to FlyDiscoPipeline() throwing an error. Another simple case is if the
    % experiment metadata indicates dead flies. That should stop the core pipeline
    % stages from running, but should not cause FlyDiscoPipeline() to throw an
    % error.
    
    fprintf('Memory usage before FlyDiscoAutomaticChecksIncoming():\n') ;
    print_matlab_memory_usage() ;
    
    %% incoming checks
    stage = 'automaticchecks_incoming';
    if is_on_or_force(doautomaticchecksincoming) ,
        forcecompute = is_force(doautomaticchecksincoming) ;
        todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticchecksincoming);
        if forcecompute || todo,
            fprintf('Running incoming automatic checks...\n');
            FlyDiscoAutomaticChecksIncoming(expdir, ...
                stage, ...
                'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
                automaticchecksincoming_params{:});
            fprintf('Memory usage after FlyDiscoAutomaticChecksIncoming():\n') ;
            print_matlab_memory_usage() ;
        end
        
        % This is just checking to make sure the 'main' file generated by
        % FlyDiscoAutomaticChecksIncoming(), usually named
        % 'automatic_checks_incoming_results.txt',
        % is present.  If it's not, it means that something has gone seriously wrong.
        [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticchecksincoming);
        if ismissingfile,
            msgs = cellfun(@(x) sprintf('Missing incoming automatic checks file %s',x),missingfiles,'UniformOutput',false) ;
            flydisco_pipeline_error(stage, msgs) ;
        end
        
        % Make sure the file usually named automatic_checks_incoming_results.txt
        % contains either "automated_pf,P" or "automated_pf,U", but not
        % "automated_pf,F" (nor anything else).
        % This ensures the the core pipeline only runs if the current run or a previous run of the ACI stage
        % indicated that it should do so.
        do_run_core_pipeline = CheckACIResultsFileContents(expdir, dataloc_params, stage) ;
        if ~do_run_core_pipeline ,
            fprintf('FlyDisco pipeline encountered issues with experiment during %s stage.  Skipping core pipeline stages.\n', stage);
        end
    else
        do_run_core_pipeline = true ;
    end
    
    
    
    if do_run_core_pipeline ,
        %% Run FlyTracker
        stage = 'flytracker' ;
        if is_on_or_force(doflytracking) ,
            forcecompute = is_force(doflytracking) ;
            todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_flytracker) ;
            if forcecompute || todo ,
                fprintf('Running FlyTracker...\n');
                FlyTrackerWrapperForFlyDisco(expdir, settingsdir, analysis_protocol, dataloc_params, forcecompute) ;
                fprintf('Memory usage after FlyTrackerWrapperForFlyDisco():\n') ;
                print_matlab_memory_usage() ;
            end
            
            % make sure flytracker files exist
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_flytracker);
            if ismissingfile,
                msgs = cellfun(@(x) sprintf('Missing FlyTracker file %s',x),missingfiles,'UniformOutput',false);
                flydisco_pipeline_error(stage, msgs) ;
            end
            
        end
        
        
        
        %% registration
        stage = 'registration';
        if is_on_or_force(doregistration) ,
            forcecompute = is_force(doregistration) ;
            todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_registration);
            if forcecompute || todo,
                fprintf('Registering tracks...\n');
                
                FlyDiscoRegisterTrx(expdir,...
                    'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
                    registration_params{:});
                fprintf('Memory usage after FlyDiscoRegisterTrx():\n') ;
                print_matlab_memory_usage() ;
            end
            
            % make sure registration files exist
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_registration);
            if ismissingfile,
                msgs = cellfun(@(x) sprintf('Missing registration file %s',x),missingfiles,'UniformOutput',false);
                flydisco_pipeline_error(stage, msgs) ;
            end
            
        end
        
        
        
        %% led indicator detection on/off
        stage = 'ledonoffdetection';
        if is_on_or_force(doledonoffdetection) ,
            forcecompute = is_force(doledonoffdetection) ;
            todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_ledonoffdetection);
            if forcecompute || todo,
                fprintf('Detecting LED on/off transitions...\n');
                FlyDiscoDectectIndicatorLedOnOff(expdir,...
                    'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
                    registration_params{:});
                fprintf('Memory usage after FlyDiscoDectectIndicatorLedOnOff():\n') ;
                print_matlab_memory_usage() ;
            end
            % make sure leddetection files exist requiredfiles_ledonoffdetection
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_ledonoffdetection);
            if ismissingfile,
                msgs = cellfun(@(x) sprintf('Missing registration file %s',x),missingfiles,'UniformOutput',false);
                flydisco_pipeline_error(stage, msgs) ;
            end
        end
        
        
        
%         %% Wing tracking and choose orientations
%         stage = 'trackwings';
%         if is_on_or_force(dotrackwings) ,
%             forcecompute = is_force(dotrackwings) ;
%             todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_wingtracking);
%             if forcecompute || todo,
%                 fprintf('Running wing tracking...\n');
%                 FlyTracker2WingTracking(expdir, ...
%                     'dataloc_params', dataloc_params, ...
%                     'settingsdir', settingsdir, ...
%                     'analysis_protocol', analysis_protocol) ;
%                 fprintf('Memory usage after FlyTracker2WingTracking():\n') ;
%                 print_matlab_memory_usage() ;
%             end
%             % make sure sexclassification files exist
%             [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_wingtracking);
%             if ismissingfile,
%                 msgs = cellfun(@(x) sprintf('Missing wing tracking file %s',x),missingfiles,'UniformOutput',false);
%                 flydisco_pipeline_error(stage, msgs) ;
%             end
%         end
        
        %% sex classification
        stage = 'sexclassification';
        if is_on_or_force(dosexclassification) ,
            forcecompute = is_force(dosexclassification) ;
            todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_sexclassification);
            if forcecompute || todo,
                fprintf('Running sex classification...\n');
                FlyDiscoClassifySex(expdir,...
                    'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
                    sexclassification_params{:});
                fprintf('Memory usage after FlyDiscoClassifySex():\n') ;
                print_matlab_memory_usage() ;
            end
            
            % make sure sexclassification files exist
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_sexclassification);
            if ismissingfile,
                msgs = cellfun(@(x) sprintf('Missing sex classification file %s',x),missingfiles,'UniformOutput',false);
                flydisco_pipeline_error(stage, msgs) ;
            end
            
        end
        
        %% compute locomotor and social per-frame features
        stage = 'computeperframefeatures';
        if is_on_or_force(docomputeperframefeatures) ,
            forcecompute = is_force(docomputeperframefeatures) ;
            
            requiredfiles_computeperframefeatures = ...
                FlyDiscoReplaceRequiredFilePlaceholders(requiredfiles_computeperframefeatures, stage, expdir, analysis_protocol_folder_path, ...
                                                        dataloc_params, computeperframefeatures_params) ;
                        
            todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframefeatures);
            if forcecompute || todo,
                fprintf('Computing per-frame features...\n');
                FlyDiscoComputePerFrameFeatures(expdir,...
                    'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
                    'forcecompute',forcecompute,...
                    computeperframefeatures_params{:});
                fprintf('Memory usage after FlyDiscoComputePerFrameFeatures():\n') ;
                print_matlab_memory_usage() ;
            end
            
            % make sure computeperframefeatures files exist
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframefeatures);
            if ismissingfile,
                msgs = cellfun(@(x) sprintf('Missing per-frame features file %s',x),missingfiles,'UniformOutput',false);
                fprintf('Computing per-frame features failed:\n');
                flydisco_pipeline_error(stage, msgs) ;
            end
            
        end
        
        %% compute perframestats
        stage = 'computeperframestats';
        if is_on_or_force(docomputeperframestats)
            forcecompute = is_force(docomputeperframestats);
            
            todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframestats);
            if forcecompute || todo,
                fprintf('ComputePerFrameStats...\n');
                stage_additional_arguments = analysis_parameters.computeperframestats_params ;
                FlyDiscoComputePerFrameStats(expdir,...
                    'settingsdir',settingsdir, ...
                    'analysis_protocol',analysis_protocol, ...
                    stage_additional_arguments{:});
                fprintf('Memory usage after FlyDiscoComputePerFrameStats():\n') ;
                print_matlab_memory_usage() ;
            end
            
            % make sure computeperframestats files exist
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframestats);
            if ismissingfile,
                msgs = cellfun(@(x) sprintf('Missing computeperframestats file %s',x),missingfiles,'UniformOutput',false);
                flydisco_pipeline_error(stage, msgs) ;
            end
            
        end
        
        
        %% compute hog hof per-frame features
        %flies_hoghof_hs_notaligned -- flow computed using Horn-Schunck. The current frame and next frame are not aligned in any way.
        
        stage = 'computehoghofperframefeatures';
        
        if is_on_or_force(docomputehoghofperframefeatures) ,
            forcecompute = is_force(docomputehoghofperframefeatures) ;
            requiredfiles_computehoghofperframefeatures = fullfile(dataloc_params.perframedir,requiredfiles_computehoghofperframefeatures);
            todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computehoghofperframefeatures);
            if forcecompute || todo,
                pwdprev = pwd;
                datalocparamsfilestr = 'dataloc_params.txt';
                datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
                dataloc_params = ReadParams(datalocparamsfile);
                trxfilestr = fullfile(expdir,dataloc_params.trxfilestr);
                movfilestr = fullfile(expdir,dataloc_params.moviefilestr);
                spacetimefeaturesdir = fileparts(which('preparePerFrameFtrs'));
                cd (spacetimefeaturesdir);
                fprintf('Computing HOG/HOF per-frame features...\n');
                preparePerFrameFtrs(movfilestr,trxfilestr,false,false);
                cd(pwdprev);
                fprintf('Memory usage after preparePerFrameFtrs():\n') ;
                print_matlab_memory_usage() ;                
            end
            % make sure computeperframefeatures files exist
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframefeatures);
            if ismissingfile,
                msgs = cellfun(@(x) sprintf('Missing HOG/HOF per-frame features file %s',x),missingfiles,'UniformOutput',false);
                flydisco_pipeline_error(stage, msgs) ;
            end
            
        end
        
        %% behavior detection
        stage = 'jaabadetect';  %#ok<NASGU>
        if is_on_or_force(dojaabadetect),
            forcecompute = is_force(dojaabadetect) ;
            JAABADetectWrapper(expdir, settingsdir, analysis_protocol, forcecompute) ;
            fprintf('Memory usage after JAABADetectWrapper():\n') ;
            print_matlab_memory_usage() ;
        end
        
        %% Make per-frame statistics plots, culminating in stats.html, which shows all of them
        stage = 'plotperframestats' ;
        stage_function_name = 'FlyDiscoPlotPerFrameStats2' ;
        FlyDiscoPipelineStage(...
            expdir, stage, doplotperframestats, dataloc_params, requiredfiles_plotperframestats, settingsdir, analysis_protocol, ...
            stage_function_name, plotperframestats_params) ;
        
        %% make results movie
        stage = 'ctraxresultsmovie';
        if is_on_or_force(domakectraxresultsmovie) ,
            forcecompute = is_force(domakectraxresultsmovie) ;
            
            requiredfiles_makectraxresultsmovie = ...
                FlyDiscoReplaceRequiredFilePlaceholders(requiredfiles_makectraxresultsmovie, stage, expdir, analysis_protocol_folder_path, ...
                                                        dataloc_params, plotperframestats_params) ;
            
            todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makectraxresultsmovie);
            if forcecompute || todo ,
                fprintf('Making Ctrax results movie...\n');
                FlyDiscoMakeCtraxResultsMovie(expdir,...
                    'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
                    makectraxresultsmovie_params{:});
                fprintf('Memory usage after FlyDiscoMakeCtraxResultsMovie():\n') ;
                print_matlab_memory_usage() ;
            end
            
            % make sure makectraxresultsmovie files exist
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makectraxresultsmovie);
            if ismissingfile,
                msgs = cellfun(@(x) sprintf('Missing Ctrax results movie file %s',x),missingfiles,'UniformOutput',false);
                flydisco_pipeline_error(stage, msgs) ;
            end
            
        end
    end  % if do_run_core_pipeline

    fprintf('Memory usage before FlyDiscoAutomaticChecksComplete():\n') ;
    print_matlab_memory_usage() ;
    
    %% automaticchecks_complete
    % This would normally be turned off for a goldblum run, and we'd run it in the
    % caboose stage.
    stage = 'automaticchecks_complete';
    if is_on_or_force(doautomaticcheckscomplete) ,
        forcecompute = is_force(doautomaticcheckscomplete) ;
        todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticcheckscomplete);
        if forcecompute || todo,
            fprintf('Running completion automatic checks...\n');
            FlyDiscoAutomaticChecksComplete(expdir,...
                'settingsdir',settingsdir, ...
                'analysis_protocol',analysis_protocol,...
                automaticcheckscomplete_params{:});
            fprintf('Memory usage after FlyDiscoAutomaticChecksComplete():\n') ;
            print_matlab_memory_usage() ;            
        end
        
        % This is just checking to make sure the 'main' file generated by
        % FlyDiscoAutomaticChecksComplete(), usually named
        % 'automatic_checks_complete_results.txt',
        % is present.  If it's not, it means that something has gone seriously wrong.
        [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticcheckscomplete);
        if ismissingfile,
            msgs = cellfun(@(x) sprintf('Missing completion automatic checks file %s',x),missingfiles,'UniformOutput',false);
            flydisco_pipeline_error(stage, msgs) ;
        end
    end
        
    %% If get here, analysis has completed successfully
    fprintf('Analysis pipeline completed!\n');    
    
    fprintf('Memory usage just before return from FlyDiscoPipeline():\n') ;
    print_matlab_memory_usage() ;        
end
