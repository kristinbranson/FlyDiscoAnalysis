function modpath()
  % Add needed libraries to Matlab path
  
  path_to_this_script = mfilename('fullpath') ;
  path_to_this_folder = fileparts(path_to_this_script) ;
  
  addpath(fullfile(path_to_this_folder, 'JAABA', 'filehandling')) ;
  addpath(fullfile(path_to_this_folder, 'JAABA', 'misc')) ;
  addpath(fullfile(path_to_this_folder, 'simplewing')) ;
  addpath(fullfile(path_to_this_folder, 'hmm')) ;
  %addpath(fullfile(path_to_parent_folder, 'flySpaceTimeFeatures')) ;
  addpath(fullfile(path_to_this_folder, 'JAABA', 'perframe')) ;
  
  % Run the FlyTracker modpath script
  flytracker_folder_path = fullfile(path_to_this_folder, 'FlyTracker') ;
  flytracker_modpath_script_path = fullfile(flytracker_folder_path, 'modpath.m') ;
  run(flytracker_modpath_script_path) ;  
  
  % Finally, add this folder itself
  addpath(path_to_this_folder) ;  % do this so that we don't have to stay in this folder
  % Also useful b/c this file may be called from elsewhere
end
