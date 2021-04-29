function FlyTrackerWrapperForFlyDisco(expdir, settingsdir, analysis_protocol, dataloc_params, do_force_tracking)
  % Process arguments
  if nargin<5 || isempty(do_force_tracking) ,
    do_force_tracking = false ;
  end    
    
  % Determine paths to needed files
  flytracker_parent_calibration_file_name = dataloc_params.flytrackerparentcalibrationstr ;
  flytracker_parent_calibration_file_path = fullfile(settingsdir, analysis_protocol, flytracker_parent_calibration_file_name) ;
  flytracker_calibration_file_name = dataloc_params.flytrackercalibrationstr ;
  flytracker_calibration_file_path = fullfile(expdir, flytracker_calibration_file_name) ;
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
  default_options.fr_sample = 200 ;
  default_options.n_flies_is_max = true;
  
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
  
  % Override certain options, since we get those directly from arguments or
  % dataloc_params
  options.f_parent_calib = flytracker_parent_calibration_file_path ;  
  options.force_tracking = do_force_tracking ;
  
  % Get the arena radius from the registration parameters, and stuff it into the
  % options. (tracker() is such that this value will override any value in the parent calibration).
  registration_parameters_file_name = dataloc_params.registrationparamsfilestr ;
  registration_parameters_file_path = fullfile(settingsdir, analysis_protocol, registration_parameters_file_name) ;
  registration_parameters = ReadParams(registration_parameters_file_path) ;
  arena_r_mm = registration_parameters.circleRadius_mm ;
  options.arena_r_mm = arena_r_mm ;
  
  % Run the tracker proper
  tracker([], options, flytracker_calibration_file_path, video_file_path) ;
end
