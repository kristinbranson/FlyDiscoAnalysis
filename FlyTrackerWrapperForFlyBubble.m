function FlyTrackerWrapperForFlyBubble(expdir, settingsdir, analysis_protocol, dataloc_params)

  video_file_name = dataloc_params.moviefilestr ;
  video_file_path = fullfile(expdir, video_file_name) ;
  
  num_cores = maxNumCompThreads;
  
  options.num_cores   = num_cores;
  options.num_chunks = options.num_cores*2;
  options.save_JAABA  = true;
  options.save_xls    = false;
  options.save_seg    = false;
  options.f_parent_calib = fullfile(settingsdir, analysis_protocol, 'calibration.mat') ;
  options.force_calib = true;
  options.expdir_naming = true;
  options.fr_sample = 200;
  
  vinfo = video_open(video_file_path);
  tracker([], options, [], vinfo);
  
  % Convert the output to ctrax form
  [~,video_base_name] = fileparts(video_file_name) ;
  flytracker_jaaba_trx_file_path = fullfile(expdir, [video_base_name '_JAABA'], 'trx.mat') ;
  ctrax_results_mat_path = fullfile(expdir, dataloc_params.ctraxfilestr) ;  
  ctrax_results_file_from_flytracker_jaaba_trx_file(ctrax_results_mat_path, flytracker_jaaba_trx_file_path)
end
