function result = determine_downstream_tracks_file_path(expdir, dataloc_params, do_run)
% Determines the path to the FlyTracker-style tracks file for stages downstream of addpflies

if is_on_or_force(do_run.addpflies) ,
  addpflies_ft_tracks_output_file_name = determine_addpflies_ft_tracks_output_file_name(dataloc_params) ;
  addpflies_ft_tracks_output_file_path = fullfile(expdir, addpflies_ft_tracks_output_file_name) ;
  result = addpflies_ft_tracks_output_file_path ;
else
  result = fullfile(expdir, dataloc_params.flytrackertrackstr) ;
end

end
