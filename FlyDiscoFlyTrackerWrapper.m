function FlyDiscoFlyTrackerWrapper(expdir, varargin)
  % Process arguments
  % NB: Leave the dataloc_params optional argument in place, b/c KB needs it.  -- ALT, 2024-08-19
  [analysis_protocol, settingsdir, forcecompute, debug, ~, dataloc_params] = ...
    myparse(varargin,...
            'analysis_protocol','current_bubble',...
            'settingsdir', default_settings_folder_path(), ...
            'forcecompute', false, ...
            'debug',false, ...
            'do_run',[],...
            'dataloc_params',[]) ;  %#ok<ASGLU> 
    
  % Read in the dataloc_params
  if isempty(dataloc_params),
    datalocparamsfile = fullfile(settingsdir, analysis_protocol, 'dataloc_params.txt') ;
    dataloc_params = ReadParams(datalocparamsfile) ;
  end

  % Determine paths to needed files
  flytracker_parent_calibration_file_name = dataloc_params.flytrackerparentcalibrationstr ;
  flytracker_parent_calibration_file_path = fullfile(settingsdir, analysis_protocol, flytracker_parent_calibration_file_name) ;
  flytracker_calibration_file_name = dataloc_params.flytrackercalibrationstr ;
  video_file_name = dataloc_params.moviefilestr ;
  video_file_path = fullfile(expdir, video_file_name) ;

  % Set default options
  default_num_cores = 2 ;  % We fix this for efficiency in production, and also reproducibility
  default_options = struct() ;
  default_options.num_cores   = default_num_cores ;
  default_options.num_chunks  = default_num_cores*2 ;
  default_options.save_JAABA  = true ;
  default_options.save_xls    = false ;
  default_options.save_seg    = false ;
  default_options.force_calib = true ;
  default_options.expdir_naming = true ;
  default_options.fr_samp = 200 ;  % Max number of frames to use when computing background model
  default_options.n_flies_is_max = true;
  default_options.do_compute_relative_features = false;  
    % Normally, when FlyTracker computes JAABA input features, it treats the case
    % of n==2 flies per chamber as special, and computes relative features like
    % 'dist_to_other' and 'angle_between' only in the n==2 case.
    % do_compute_relative_features==true will yield this default behavior.
    % do_compute_relative_features==false implies that the case
    % of n==2 flies per chamber will *not* cause relative features to be computed.
    %  I.e. if this is false, n==2 is treated like any other n.
  default_options.save_JAABA_trk_mat_file = true ;  % whether to save the JABBA/CTRAX trk.mat.
  default_options.save_JAABA_perframe_features = false ;  % whether to save the JABBA perframe feature .mat files.
  default_options.save_JAABA_movie_link_or_copy = true ;  % whether to create a link/copy of the movie in the movie_JAABA folder.
  
  % Read the options file, if dataloc param specifies it, and it exists
  if isfield(dataloc_params, 'flytrackeroptionsstr') ,
    flytracker_options_file_name = dataloc_params.flytrackeroptionsstr ;
    flytracker_options_file_path = fullfile(settingsdir, analysis_protocol, flytracker_options_file_name) ;    
    if exist(flytracker_options_file_path, 'file') ,
      options_from_file = ReadParams(flytracker_options_file_path) ;
    else
      options_from_file = struct() ;
    end
  else
      options_from_file = struct() ;
  end
  
  % Merge options
  options = merge_structs(default_options, options_from_file) ;
%   % A hack for running locally in some situations
%   options.num_cores = 8 ;
%   options.num_chunks = 16 ;
%   options.force_calib = false ;
  
  % Override certain options, since we get those directly from arguments or
  % dataloc_params
  %options.f_parent_calib = flytracker_parent_calibration_file_path ;  
  options.do_recompute_tracking = forcecompute ;
  
  % Get the arena radius from the registration parameters, and stuff it into the
  % options. (tracker() is such that this value will override any value in the parent calibration).
  % The "if isfield(dataloc_params,'registrationparamsfilestr')" check is needed
  % to support non-circular arenas used by external code.  -- ALT, 2024-08-18
  if isfield(dataloc_params,'registrationparamsfilestr'),
    registration_parameters_file_name = dataloc_params.registrationparamsfilestr ;
    registration_parameters_file_path = fullfile(settingsdir, analysis_protocol, registration_parameters_file_name) ;
    registration_parameters = ReadParams(registration_parameters_file_path) ;
    arena_r_mm = registration_parameters.circleRadius_mm ;
    options.arena_r_mm = arena_r_mm ;
  end
  
  % Run the tracker proper
  fda_batch_track_single_video(expdir, video_file_path, flytracker_parent_calibration_file_path, options, flytracker_calibration_file_name)  
end
