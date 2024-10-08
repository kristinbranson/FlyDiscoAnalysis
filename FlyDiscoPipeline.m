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
    %   points to settings-internal/, which *is* included in the source
    %   repository.
    %
    %   FlyDiscoPipeline(expdir, 'settingsdir', settingsdir) looks in the given
    %   settingsdir for analysis-protocol folders, rather than the default settings
    %   folder described above. 
    %
    %   FlyDiscoPipeline(expdir, 'analysis_protocol', analysis_protocol) uses the
    %   analysis-protocol folder specified instead of the analysis-protocol folder
    %   indicated by the experiment metadata.  analysis_protocol should be a string
    %   that corresponds to a subfolder name within the settings folder.
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
    %       addpflies
    %           Adds 'projector flies' to the tracks.
    %           Default: 'off'. 
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
    %       locomotionmetrics
    %           Compute locomotion metrics.
    %           Default: 'off'. 
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
    %   All stages support similar lists of additional parameters.  For
    %   each such stage, the name of the relevant parameter to FlyDiscoPipeline() is
    %   the stage name followed by '_params'.  See the documentation for the
    %   function that implements each stage for details about what parameters are
    %   supported.  Some of the stages, and the function
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
    %   In fact, all the optional parameters described above (all the
    %   do<stage_name> and <stage_name>_params parameters) are determined by
    %   starting with a set of default parameter values, overriding these with
    %   values set in the analysis-protocol folder, and finally overriding these
    %   with values provided in the argument list to FlyDiscoPipeline().
    %
    %   The default value for all these parameters is specified by the function
    %   FlyDiscoPipelineDefaultAnalysisParameters().  The values returned by that
    %   function are (possibly) overridden by values specified by the file
    %   analysis-protocol-parameters.txt in the analysis-protocol folder.  That file
    %   specifies a set of name-value pairs, and is read using the function
    %   ReadParams().  The values specified in that file override those specified in
    %   FlyDiscoPipelineDefaultAnalysisParameters().  Finally, any of these values can
    %   in turn be overridden by providing a name-value pair in the argument list to
    %   FlyDiscoPipeline().

    
    % Report the Matlab version
    matlab_ver_string = version() ;
    fprintf('Matlab version:\n%s\n\n', matlab_ver_string) ;

    % Get info about the state of the repo, output to log
    fprintf('FlyDiscoAnalysis repository state:\n')
    this_script_path = mfilename('fullpath') ;
    source_folder_path = fileparts(this_script_path) ;
    git_report = get_git_report(source_folder_path) ;
    fprintf('%s', git_report) ;

    % Process varargin parameters---If there's one optional arg and it's a struct,
    % then each field is an optional arg name, and each field's value is, well, a
    % value.  Otherwise, treat the optional arg as name-value pairs and collect
    % them up into a struct.
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
    
    % Read in the analysis parameters, dataloc_params
    [analysis_parameters_according_to_protocol, dataloc_params, analysis_protocol_folder_path] = ...
        readIntermediateAnalysisParameters(settingsdir, analysis_protocol) ;     % analysis params according to defaults and the analysis-protocol folder

    % Print the canonical path to the analysis folder
    canonical_analysis_protocol_folder_path = realpath(absolute_filename(analysis_protocol_folder_path)) ;
    fprintf('Canonical path to analysis protocol folder is:\n  %s\n\n', canonical_analysis_protocol_folder_path) ;
    
    % Print the canonical path to the experiment folder
    canonical_experiment_folder_path = realpath(absolute_filename(expdir)) ;
    fprintf('Canonical path to experiment folder is:\n  %s\n\n', canonical_experiment_folder_path) ;
        
    % Merge in the argument parameters
    analysis_parameters = merge_structs(analysis_parameters_according_to_protocol, argument_parameters) ;
    do_run = do_run_from_analysis_parameters(analysis_parameters) ;  % whether or not to run each stage
    
    % Compute the required output files for each stage
    required_file_names_from_stage_name = FlyDiscoComputeRequiredFileNamesForAllStages(analysis_protocol_folder_path, dataloc_params, expdir, do_run) ;

    % Assign a few params to variables
    cluster_billing_account_name = lookup_in_struct(analysis_parameters, 'cluster_billing_account_name') ;
    sshhost = lookup_in_struct(analysis_parameters, 'sshhost') ;
    debug = lookup_in_struct(analysis_parameters, 'debug') ;
    do_try = lookup_in_struct(analysis_parameters, 'do_try', true) ;
      % Whether or not we wrap the running of each stage in a try-catch clause.
      % Should be true in production runs to ensure that the final ACC stage runs
      % even if a core stage throws.

    % Handle the per-stage parameters
    stage_function_from_stage_name = FlyDiscoStageFunctionFromName() ;
    stage_name_from_stage_index = fieldnames(stage_function_from_stage_name) ;
    stage_count = numel(stage_name_from_stage_index) ;
    stage_additional_arguments_from_stage_name = struct() ;
    for stage_index = 1 : stage_count ,
        stage_name = stage_name_from_stage_index{stage_index} ;
        analysis_parameters_field_name = horzcat(stage_name, '_params') ;
        stage_additional_arguments = lookup_in_struct(analysis_parameters, analysis_parameters_field_name, cell(1,0)) ;
        stage_additional_arguments_from_stage_name.(stage_name) = stage_additional_arguments ;
    end
    
    % Need to supply custom arguments for a few stages
    stage_additional_arguments_from_stage_name.apt = ...
      horzcat(stage_additional_arguments_from_stage_name.apt, ...
              {'cluster_billing_account_name', cluster_billing_account_name, ...
               'sshhost', sshhost}) ;
    stage_additional_arguments_from_stage_name.automaticcheckscomplete = ...
      horzcat(stage_additional_arguments_from_stage_name.automaticcheckscomplete, ...
              {'required_file_names_from_stage_name', required_file_names_from_stage_name}) ;

    % Print the settings values in use
    variable_names_to_print = ...
        { 'settingsdir', ...
          'analysis_protocol'}' ;
    fprintf('Settings values in FlyDiscoPipeline():\n') ;
    for i = 1 : length(variable_names_to_print) ,
        variable_name = variable_names_to_print{i} ;
        value = eval(variable_name) ;
        fprintf('  %s: %s\n', variable_name, char_array_from_value(value)) ;
    end
    for stage_index = 1 : numel(stage_name_from_stage_index) ,
      stage_name = stage_name_from_stage_index{stage_index} ;
      value = do_run.(stage_name) ;
      fprintf('  do_run.%s: %s\n', stage_name, value) ;      
    end
    fprintf('\n') ;
    
    %% check that experiment exists
    if ~exist(expdir, 'dir') ,
        error('Experiment directory %s does not exist', expdir) ;
    end
    
    
    
    %% Run the stages
    FlyDiscoPipelineCore(expdir, ...
                         stage_name_from_stage_index, ...
                         do_run, ...
                         required_file_names_from_stage_name, ...
                         debug, ...
                         settingsdir, ...
                         analysis_protocol, ...
                         stage_function_from_stage_name, ...
                         stage_additional_arguments_from_stage_name, ...
                         dataloc_params, ...
                         do_try) ;
end  % function
