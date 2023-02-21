function FlyDiscoPipeline(expdir, varargin)
    %FlyDiscoPipeline  Runs the FlyDisco pipeline.
    %   FlyDiscoPipeline(expdir) runs the FlyDisco pipeline on the experiment folder
    %   expdir.  The "analysis protocol" to be used, including what FlyDisco stages to
    %   run, is determined from the screen_type field in the the experiment
    %   metadata.  For example, one possible screen_type is 'FlyBowlRGBbasic'.
    %   The experiment metadata is stored in the file metadata.xml (case
    %   insensitive) within the experiment folder.
    %
    %   An analysis protocol is specified by the files within a folder called an
    %   "analysis-protocol folder".  These files specify parameters to be used for
    %   each stage of the analysis, such as what the minimum and maximum number of
    %   flies can be, and what the shape of the arena is.  (A number of
    %   analysis-protocol folders are provided in the software repository that
    %   contains the source code for FlyDiscoPipeline(), within a folder named
    %   "settings-internal" at the top level of the repository.)  The screen_type
    %   of the experiment is used to determine which analysis-protocol folder is
    %   used to analyze that experiment.  First, a folder named
    %   'current_'+screen_type is checked for (where "+" here represents
    %   concatenation).  If present, that folder is used as the analysis-protocol
    %   folder.  If absent, a folder with the same name as the screen_type is
    %   checked for, and used if present.  For instance, one analysis-protocol
    %   folder is named 'current_FlyBowlRGBbasic'. (Which as of this writing is
    %   actually a symbolic link to an analysis-protocol folder named
    %   '20210913_FlyBowlRGBbasic'.)
    %
    %   The expdir is the only required input to FlyDiscoPipeline.  Addtional
    %   optional arguments are supported as name-value pairs.  Any combination of
    %   these name-value pairs is generally permitted, and their order does not
    %   matter.
    %
    %   The folder that is checked for analysis-protocol folders is called a
    %   "settings" folder.  The default settings folder is named "settings"
    %   (predictably), and is a sibling to the FlyDiscoPipeline.m source code file.
    %   Note, however, that no "settings/" folder is included in the source
    %   respository.  This is to prevent unintended usage of the wrong settings/
    %   folder.  If convenient, a softlink named "settings" can be created that
    %   points to settings-for-testing/, which *is* included in the source
    %   repository.
    %
    %   FlyDiscoPipeline(expdir, 'settingsdir', settingsdir) looks in the given
    %   settingsdir for analysis-protocol folders, rather than the default settings
    %   folder described above. 
    %
    %   FlyDiscoPipeline(expdir, 'analysis_protocol', analysis_protocol) uses the
    %   analysis-protocol folder specified instead of the analysis-protocol folder
    %   indicated by the experiment metadata.  analysis_protocol should be a string
    %   that corresponds to a subfolder name within the settings/ folder.
    %
    %   The operation of FlyDiscoPipeline() consists of a number of stages.  Each of
    %   these stages can be turned on or off independently.  (Although some stages
    %   will fail if earlier stages that feed into them are turned off.)  The stages
    %   are:
    %
    %       automaticchecksincoming
    %           Makes sure all the raw experiment files are present. 
    %           Default: 'force'. 
    %       flytracking
    %           Tracks the flies in the video, computing x,y coordinates and heading
    %           angle for each in each frame.
    %           Default: 'on'. 
    %       registration
    %           Determines where the boundaries of the chamber are, and where each
    %           fly is relative to those boundaries.  Also converts fly positions in
    %           pixels to positions in millimeters, relative to the center of the
    %           chamber.
    %           Default: 'on'. 
    %       ledonoffdetection
    %           Determines when each LED turns on and off.
    %           Default: 'on'. 
    %       sexclassification
    %           Determines the sex of each fly.
    %           Default: 'on'. 
    %       computeperframefeatures
    %           Computes a set of per-frame features for each fly, specified by the
    %           analysis-protocol file.  These may include things like the fly's
    %           velocity and acceleration in each frame.
    %           Default: 'on'. 
    %       computeperframestats
    %           Compute statistics of the per-frame features across all frames of
    %           the video.
    %           Default: 'off'. 
    %       computehoghofperframefeatures
    %           Compute Histogram-of-Oriented-Gradient (HOG) features for each
    %           frame, for each fly.  Also HOF features. Not tested, makes large files
    %           Default: 'off'. 
    %       jaabadetect
    %           Use JAABA to classify the behavior of each fly in each frame
    %           according to a set of classes.
    %           Default: 'on'. 
    %       plotperframestats
    %           Make histograms and other plots of the per-frame features across
    %           frames.
    %           Default: 'off'. 
    %       makectraxresultsmovie
    %           Make a .mp4 movie summarizing the results of the tracking.  (The
    %           CTRAX in the name is vestigial---CTRAX is not used anymore.)
    %           Default: 'on'. 
    %       apt
    %           Run APT, the Advanced Part Tracker, on the video.
    %           Default: 'off'. 
    %       makeaptresultmovie
    %           Make a .mp4 movie summarizing the results of APT pose tracking. 
    %            Default: 'off'.
    %       automaticcheckscomplete
    %           Perform a set of final checks that the analysis of the experiment
    %           has been completed successfully.
    %           Default: 'force'. 
    %   Default values are specified in FlyDiscoPipelineDefaultAnalysisParameters.m. 
    %
    %   FlyDiscoPipeline(expdir, 'doautomaticchecksincoming', offonforce), where
    %   offonforce is either 'off', 'on', or 'force', specifies whether or not to run the
    %   automaticchecksincoming stage of the pipeline.  'off' means do not run this
    %   stage, 'on' means run it if the stage's output files are missing from the
    %   experiment folder, and 'force' means run it regardless of whether the output
    %   files are present.
    %
    %   All other stages can be turned on/off (or forced) in a similar fashion.  The
    %   name of the relevant parameter is always 'do' followed by the stage name.
    %
    %   (All text below this point documents rarely-used features, and can be safely
    %   ignored unless you are an advanced user.)
    %
    %   FlyDiscoPipeline(expdir, 'automaticchecksincoming_params', params) passes
    %   the given params (typically a cell array of name-value pairs) to the
    %   automaticchecksincoming stage, which is implemented by the function
    %   FlyDiscoAutomaticChecksIncoming().  See the documentation for that function
    %   for details about what parameters are supported.
    %
    %   Several other stages support similar lists of additional parameters.  For
    %   each such stage, the name of the relevant parameter to FlyDiscoPipeine() is
    %   the stage name followed by '_params'.  See the documentation for the
    %   function that implements each stage for details about what parameters are
    %   supported.  The stages that support a '_params' parameter, and the function
    %   implementing the respective stage, are:
    %
    %       automaticchecksincoming: FlyDiscoAutomaticChecksIncoming()
    %       sexclassification: FlyDiscoClassifySex()
    %       computeperframefeatures: FlyDiscoComputePerFrameFeatures()
    %       computeperframestats: FlyDiscoComputePerFrameStats()
    %       plotperframestats: FlyDiscoPlotPerFrameStats()
    %       makectraxresultsmovie: FlyDiscoMakeCtraxResultsMovie()
    %       apt: FlyDiscoAPTTrack()
    %       makeaptresultmovie: FlyDiscoMakeAPTResultsMovie()
    %       automaticcheckscomplete: FlyDiscoAutomaticChecksComplete()
    %
    %   In fact, all the optional parameters described above (besides 'settingsdir'
    %   and 'analysis_protocol') are determined by starting with a set of default
    %   parameter values, overriding these with values set in the analysis-protocol
    %   folder, and finally overriding these with values provided in the argument
    %   list to FlyDiscoPipeline().  The full list of parameters that use this
    %   three-level scheme is:
    %
    %       datalocparamsfilestr
    %       automaticchecksincoming_params
    %       registration_params
    %       sexclassification_params
    %       computeperframefeatures_params
    %       computeperframestats_params
    %       plotperframestats_params
    %       makectraxresultsmovie_params
    %       aptresultsmovie_params
    %       automaticcheckscomplete_params
    %       doautomaticchecksincoming
    %       doflytracking
    %       doregistration
    %       doledonoffdetection
    %       dosexclassification
    %       docomputeperframefeatures
    %       docomputehoghofperframefeatures
    %       dojaabadetect
    %       docomputeperframestats
    %       doplotperframestats
    %       domakectraxresultsmovie
    %       doapt
    %       domakeaptresultsmovie
    %       doautomaticcheckscomplete
    %       requiredfiles_automaticchecksincoming
    %       requiredfiles_flytracker
    %       requiredfiles_registration
    %       requiredfiles_ledonoffdetection
    %       requiredfiles_sexclassification
    %       requiredfiles_computeperframefeatures
    %       requiredfiles_computeperframestats
    %       requiredfiles_computehoghofperframefeatures
    %       requiredfiles_plotperframestats
    %       requiredfiles_makectraxresultsmovie
    %       requiredfiles_apt    
    %       requiredfiles_makeaptresultsmovie
    %       requiredfiles_automaticcheckscomplete
    %
    % The default value for all these parameters is specified by the function
    % FlyDiscoPipelineDefaultAnalysisParameters().  The values returned by that
    % function are (possibly) overridden by values specified by the file
    % analysis-protocol-parameters.txt in the analysis-protocol folder.  That file
    % specifies a set of name-value pairs, and is read using the function
    % ReadParams().  The values specified in that file override those specified in
    % FlyDiscoPipelineDefaultAnalysisParameters().  Finally, any of these values can
    % in turn be overridden by providing a name-value pair in the argument list to
    % FlyDiscoPipeline().
    %
    % It should be noted that not all parameters of all stages can be set in this
    % way.  Some parameters are determined by the files in the analysis-protocol
    % folder, and cannot easily be changed without modifying the analysis-protcol
    % folder.

    
    % Report the Matlab version
    matlab_ver_string = version() ;
    fprintf('Matlab version:\n%s\n\n', matlab_ver_string) ;

    % Get info about the state of the repo, output to log
    fprintf('FlyDiscoAnalysis repository state:\n')
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
    settingsdir = lookup_in_struct(argument_parameters, 'settingsdir', default_settings_folder_path()) ;

    % If the settings folder is not part of the FDA repo, print a report about it
    internal_settings_folder_path = fullfile(source_folder_path, 'settings-internal') ;
    canonical_settings_folder_path = realpath(settingsdir) ;
    fprintf('Settings folder repository state:\n')
    if strcmp(canonical_settings_folder_path, internal_settings_folder_path) ,
        fprintf('(Using internal settings folder)\n\n\n')
    else
        settings_git_report = get_git_report(settingsdir) ;
        fprintf('%s', settings_git_report) ;        
    end

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
    computeperframestats_params = lookup_in_struct(analysis_parameters, 'computeperframestats_params') ;    
    plotperframestats_params = lookup_in_struct(analysis_parameters, 'plotperframestats_params') ;    
    makectraxresultsmovie_params = lookup_in_struct(analysis_parameters, 'makectraxresultsmovie_params') ;
    aptresultsmovie_params = lookup_in_struct(analysis_parameters, 'aptresultsmovie_params') ;
    automaticcheckscomplete_params = lookup_in_struct(analysis_parameters, 'automaticcheckscomplete_params') ;
    doautomaticchecksincoming = lookup_in_struct(analysis_parameters, 'doautomaticchecksincoming') ;
    doflytracking = lookup_in_struct(analysis_parameters, 'doflytracking') ;
    doregistration = lookup_in_struct(analysis_parameters, 'doregistration') ;
    doledonoffdetection = lookup_in_struct(analysis_parameters, 'doledonoffdetection') ;
    dosexclassification = lookup_in_struct(analysis_parameters, 'dosexclassification') ;
    docomputeperframefeatures = lookup_in_struct(analysis_parameters, 'docomputeperframefeatures') ;
    docomputehoghofperframefeatures = lookup_in_struct(analysis_parameters, 'docomputehoghofperframefeatures') ;
    dojaabadetect = lookup_in_struct(analysis_parameters, 'dojaabadetect') ;
    docomputeperframestats = lookup_in_struct(analysis_parameters, 'docomputeperframestats') ;
    doplotperframestats = lookup_in_struct(analysis_parameters, 'doplotperframestats') ;
    domakectraxresultsmovie = lookup_in_struct(analysis_parameters, 'domakectraxresultsmovie') ;
    doapt = lookup_in_struct(analysis_parameters, 'doapt') ;
    domakeaptresultsmovie = lookup_in_struct(analysis_parameters, 'domakeaptresultsmovie') ;
    doautomaticcheckscomplete = lookup_in_struct(analysis_parameters, 'doautomaticcheckscomplete') ;
    requiredfiles_automaticchecksincoming = lookup_in_struct(analysis_parameters, 'requiredfiles_automaticchecksincoming') ;
    requiredfiles_flytracker = lookup_in_struct(analysis_parameters, 'requiredfiles_flytracker') ;
    requiredfiles_registration = lookup_in_struct(analysis_parameters, 'requiredfiles_registration') ;
    requiredfiles_ledonoffdetection = lookup_in_struct(analysis_parameters, 'requiredfiles_ledonoffdetection') ;
    requiredfiles_sexclassification = lookup_in_struct(analysis_parameters, 'requiredfiles_sexclassification') ;
    requiredfiles_computeperframefeatures = lookup_in_struct(analysis_parameters, 'requiredfiles_computeperframefeatures') ;
    requiredfiles_computeperframestats = lookup_in_struct(analysis_parameters, 'requiredfiles_computeperframestats') ;
    requiredfiles_computehoghofperframefeatures = lookup_in_struct(analysis_parameters, 'requiredfiles_computehoghofperframefeatures') ;
    requiredfiles_plotperframestats = lookup_in_struct(analysis_parameters, 'requiredfiles_plotperframestats') ;    
    requiredfiles_makectraxresultsmovie = lookup_in_struct(analysis_parameters, 'requiredfiles_makectraxresultsmovie') ;
    requiredfiles_apt = lookup_in_struct(analysis_parameters, 'requiredfiles_apt') ;
    requiredfiles_makeaptresultsmovie = lookup_in_struct(analysis_parameters, 'requiredfiles_makeaptresultsmovie') ;    
    requiredfiles_automaticcheckscomplete = lookup_in_struct(analysis_parameters, 'requiredfiles_automaticcheckscomplete') ;
    cluster_billing_account_name = lookup_in_struct(analysis_parameters, 'cluster_billing_account_name') ;

    % Read in the dataloc params
    datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
    dataloc_params = ReadParams(datalocparamsfile);
    
    % Coerce all the do* variables to on/off/force
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
    doapt = coerce_to_on_off_force(doapt) ;
    domakeaptresultsmovie = coerce_to_on_off_force(domakeaptresultsmovie) ;

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
          'docomputeperframestats', ...
          'docomputehoghofperframefeatures', ...
          'dojaabadetect', ...
          'doplotperframestats', ...
          'domakectraxresultsmovie', ...
          'doapt', ...
          'domakeaptresultsmovie', ...
          'doautomaticcheckscomplete' }' ;
    fprintf('Settings values in FlyDiscoPipeline():\n') ;
    for i = 1 : length(variable_names_to_print) ,
        variable_name = variable_names_to_print{i} ;
        value = eval(variable_name) ;
        fprintf('  %s: %s\n', variable_name, char_array_from_value(value)) ;
    end
    fprintf('\n') ;
    
    % Print the canonical path to the analysis folder
    canonical_analysis_protocol_folder_path = realpath(absolute_filename(analysis_protocol_folder_path)) ;
    fprintf('Canonical path to analysis protocol folder is:\n  %s\n\n', canonical_analysis_protocol_folder_path) ;
    
    % Print the canonical path to the experiment folder
    canonical_experiment_folder_path = realpath(absolute_filename(expdir)) ;
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
    
%     fprintf('Memory usage before FlyDiscoAutomaticChecksIncoming():\n') ;
%     print_matlab_memory_usage() ;
    
    %% incoming checks
    stage = 'automaticchecksincoming';
    if is_on_or_force(doautomaticchecksincoming) ,
        forcecompute = is_force(doautomaticchecksincoming) ;
        todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticchecksincoming);
        if forcecompute || todo,
            fprintf('Running incoming automatic checks...\n');
            FlyDiscoAutomaticChecksIncoming(expdir, ...
                stage, ...
                'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
                automaticchecksincoming_params{:});
%             fprintf('Memory usage after FlyDiscoAutomaticChecksIncoming():\n') ;
%             print_matlab_memory_usage() ;
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
        stage = 'flytracking' ;
        if is_on_or_force(doflytracking) ,
            forcecompute = is_force(doflytracking) ;
            todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_flytracker) ;
            if forcecompute || todo ,
                fprintf('Running FlyTracker...\n');
                FlyTrackerWrapperForFlyDisco(expdir, settingsdir, analysis_protocol, dataloc_params, forcecompute) ;
%                 fprintf('Memory usage after FlyTrackerWrapperForFlyDisco():\n') ;
%                 print_matlab_memory_usage() ;
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
%                 fprintf('Memory usage after FlyDiscoRegisterTrx():\n') ;
%                 print_matlab_memory_usage() ;
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
%                 fprintf('Memory usage after FlyDiscoDectectIndicatorLedOnOff():\n') ;
%                 print_matlab_memory_usage() ;
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
%                 fprintf('Memory usage after FlyDiscoClassifySex():\n') ;
%                 print_matlab_memory_usage() ;
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
%                 fprintf('Memory usage after FlyDiscoComputePerFrameFeatures():\n') ;
%                 print_matlab_memory_usage() ;
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
                stage_additional_arguments = computeperframestats_params ;
                FlyDiscoComputePerFrameStats(expdir,...
                    'settingsdir',settingsdir, ...
                    'analysis_protocol',analysis_protocol, ...
                    stage_additional_arguments{:});
%                 fprintf('Memory usage after FlyDiscoComputePerFrameStats():\n') ;
%                 print_matlab_memory_usage() ;
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
                %datalocparamsfilestr = 'dataloc_params.txt';
                %datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
                %dataloc_params = ReadParams(datalocparamsfile);
                trxfilestr = fullfile(expdir,dataloc_params.trxfilestr);
                movfilestr = fullfile(expdir,dataloc_params.moviefilestr);
                spacetimefeaturesdir = fileparts(which('preparePerFrameFtrs'));
                cd (spacetimefeaturesdir);
                fprintf('Computing HOG/HOF per-frame features...\n');
                preparePerFrameFtrs(movfilestr,trxfilestr,false,false);
                cd(pwdprev);
%                 fprintf('Memory usage after preparePerFrameFtrs():\n') ;
%                 print_matlab_memory_usage() ;                
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
%             fprintf('Memory usage after JAABADetectWrapper():\n') ;
%             print_matlab_memory_usage() ;
        end
        
        %% Make per-frame statistics plots, culminating in stats.html, which shows all of them
        stage = 'plotperframestats' ;
        stage_function_name = 'FlyDiscoPlotPerFrameStats' ;
        FlyDiscoPipelineStage(...
            expdir, stage, doplotperframestats, dataloc_params, requiredfiles_plotperframestats, settingsdir, analysis_protocol, ...
            stage_function_name, plotperframestats_params) ;
        
        %% make results movie
        stage = 'makectraxresultsmovie';
        if is_on_or_force(domakectraxresultsmovie) ,
            forcecompute = is_force(domakectraxresultsmovie) ;
            
            requiredfiles_makectraxresultsmovie = ...
                FlyDiscoReplaceRequiredFilePlaceholders(requiredfiles_makectraxresultsmovie, stage, expdir, analysis_protocol_folder_path, ...
                                                        dataloc_params, makectraxresultsmovie_params) ;
            
            todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makectraxresultsmovie);
            if forcecompute || todo ,
                fprintf('Making Ctrax results movie...\n');
                FlyDiscoMakeCtraxResultsMovie(expdir, ...
                                              'settingsdir',settingsdir, ...
                                              'analysis_protocol',analysis_protocol, ...
                                              makectraxresultsmovie_params{:}) ;
            end
            
            % make sure makectraxresultsmovie files exist
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makectraxresultsmovie);
            if ismissingfile,
                msgs = cellfun(@(x) sprintf('Missing Ctrax results movie file %s',x),missingfiles,'UniformOutput',false);
                flydisco_pipeline_error(stage, msgs) ;
            end
            
        end
        
        %% Run APT
        stage = 'apt' ;
        if is_on_or_force(doapt) ,
            forcecompute = is_force(doapt) ;
            todo = CheckForMissingFiles(expdir, dataloc_params, requiredfiles_apt) ;
            if forcecompute || todo ,
                fprintf('Running APT...\n');                
                aptrepopath = fullfile(source_folder_path, 'APT') ;  % use the subrepo APT
                % Are we running on a submit host?  If not, use submit.int.janelia.org.
                [return_code, ~] = system('/usr/bin/which bsub') ;
                if return_code == 0 ,
                    submit_host_name = '' ;
                else
                    submit_host_name = 'submit.int.janelia.org' ;
                end
                [success, msgs] = ...
                    FlyDiscoAPTTrack(expdir, ...
                                     'settingsdir',settingsdir,...
                                     'analysis_protocol',analysis_protocol,...
                                     'dataloc_params', dataloc_params, ...
                                     'dryrun',false,...
                                     'verbose',1,...
                                     'dowait',true,...
                                     'waitchecktime',300,...
                                     'dooverwrite',forcecompute, ...
                                     'cluster_billing_account_name', cluster_billing_account_name, ...
                                     'aptrepopath', aptrepopath, ...
                                     'docomputemd5s', true, ...
                                     'sshhost', submit_host_name) ;
                if ~success ,
                    flydisco_pipeline_error(stage, msgs) ;
                end
            end
            
            % make sure flytracker files exist
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir, dataloc_params, requiredfiles_apt) ;
            if ismissingfile ,
                msgs = cellfun(@(x) sprintf('Missing APT file %s',x),missingfiles,'UniformOutput',false);
                flydisco_pipeline_error(stage, msgs) ;
            end
        end  % APT stage
        %% make apt results movie
        stage = 'makeaptresultsmovie';
        if is_on_or_force(domakeaptresultsmovie)
            forcecompute = is_force(domakeaptresultsmovie);
            requiredfiles_makeaptresultsmovie = ...
                FlyDiscoReplaceRequiredFilePlaceholders(requiredfiles_makeaptresultsmovie, stage, expdir, analysis_protocol_folder_path, ...
                dataloc_params, aptresultsmovie_params) ;
            todo = CheckForMissingFiles(expdir, dataloc_params, requiredfiles_makeaptresultsmovie) ;
            if forcecompute || todo
                fprintf('Making APT results movie... \n');
                FlyDiscoMakeAPTResultsMovie(expdir, ...
                    'settingsdir',settingsdir, ...
                    'analysis_protocol',analysis_protocol, ...
                    aptresultsmovie_params{:}) ;
                
            end
            % make sure makeaptresultsmovie files exist
            [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makeaptresultsmovie);
            if ismissingfile,
                msgs = cellfun(@(x) sprintf('Missing APT results movie file %s',x),missingfiles,'UniformOutput',false);
                flydisco_pipeline_error(stage, msgs) ;
            end
        end % make apt results movie
    end  % if do_run_core_pipeline

%     fprintf('Memory usage before FlyDiscoAutomaticChecksComplete():\n') ;
%     print_matlab_memory_usage() ;
    
    %% automaticchecks_complete
    % This would normally be turned off for a goldblum run, and we'd run it in the
    % caboose stage.
    stage = 'automaticcheckscomplete';
    if is_on_or_force(doautomaticcheckscomplete) ,
        forcecompute = is_force(doautomaticcheckscomplete) ;
        todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticcheckscomplete);
        if forcecompute || todo,
            fprintf('Running completion automatic checks...\n');
            FlyDiscoAutomaticChecksComplete(expdir,...
                'settingsdir',settingsdir, ...
                'analysis_protocol',analysis_protocol,...
                automaticcheckscomplete_params{:});
%             fprintf('Memory usage after FlyDiscoAutomaticChecksComplete():\n') ;
%             print_matlab_memory_usage() ;            
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
    
%     fprintf('Memory usage just before return from FlyDiscoPipeline():\n') ;
%     print_matlab_memory_usage() ;        
end
