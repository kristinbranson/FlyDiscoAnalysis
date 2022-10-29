function modpath()
  % Add needed libraries to Matlab path

  % Sort out where FlyBowlAnalysis is in the filesytem
  path_to_this_script = mfilename('fullpath') ;
  path_to_this_folder = fileparts(path_to_this_script) ;
  
  % Run the goldblum modpath
  goldblum_folder_path = fullfile(path_to_this_folder, 'goldblum') ;
  goldblum_modpath_script_path = fullfile(goldblum_folder_path, 'modpath.m') ;
  run(goldblum_modpath_script_path) ;    
  
  % Run the FlyTracker modpath script
  flytracker_folder_path = fullfile(path_to_this_folder, 'FlyTracker') ;
  flytracker_modpath_script_path = fullfile(flytracker_folder_path, 'modpath.m') ;
  run(flytracker_modpath_script_path) ;  

  % Add the JAABA stuff to the path
  jaaba_perframe_folder_path = fullfile(path_to_this_folder, 'JAABA', 'perframe') ;
  jaaba_modpath_script_path = fullfile(jaaba_perframe_folder_path, 'SetUpJAABAPath.m') ;
  run(jaaba_modpath_script_path) ;    
  
  % Add FlyDiscoAnalysis subfolders that are not their own projects
  addpath(fullfile(path_to_this_folder, 'simplewing')) ;
  addpath(fullfile(path_to_this_folder, 'hmm')) ;  
  addpath(fullfile(path_to_this_folder, 'filehandling')) ;  
  addpath(fullfile(path_to_this_folder, 'perframe')) ;

  % Add stuff intended to shadow JAABA, FlyTracker versions
  addpath(fullfile(path_to_this_folder, 'shadow')) ;
  
  % Finally, add this folder itself, so we don't have to stay in this folder
  addpath(path_to_this_folder) ;
  
  % Run code to set the parpool location appropriately
  set_parpool_job_storage_location()
end
