function modpath()
  % Add needed libraries to Matlab path
  
  path_to_this_script = mfilename('fullpath') ;
  path_to_this_folder = fileparts(path_to_this_script) ;
  path_to_parent_folder = fileparts(path_to_this_folder) ;
  addpath(fullfile(path_to_parent_folder, 'JAABA', 'filehandling')) ;
  addpath(fullfile(path_to_parent_folder, 'JAABA', 'misc')) ;
  addpath(fullfile(path_to_this_folder, 'simplewing')) ;
  addpath(fullfile(path_to_this_folder, 'hmm')) ;
  %addpath(fullfile(path_to_parent_folder, 'flySpaceTimeFeatures')) ;
  addpath(fullfile(path_to_parent_folder, 'JAABA', 'perframe')) ;
  
  % Run the FlyTracker modpath script
  % We assume FlyTracker is in a sibling folder named "FlyTracker"
  % We'll include it as a submodule once things get sorted with that.
  flytracker_folder_path = fullfile(path_to_parent_folder, 'FlyTracker') ;
  flytracker_modpath_script_path = fullfile(flytracker_folder_path, 'modpath.m') ;
  run(flytracker_modpath_script_path) ;  
end
