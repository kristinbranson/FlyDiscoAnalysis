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

  % Don't SetUpJAABAPath(), that will shadow the FlyBowlAnalysis compute_*
  % functions.
  % We'll add that stuff to the path in DetectJAABA(), and tear it down just
  % before exit
  jaaba_perframe_folder_path = fullfile(path_to_this_folder, 'JAABA', 'perframe') ;
  addpath(jaaba_perframe_folder_path) ;
  addpath(fullfile(path_to_this_folder, 'JAABA', 'filehandling')) ;  % ok, add this too b/c need load_tracks(), at least
  addpath(fullfile(path_to_this_folder, 'JAABA', 'misc')) ;  % ok, add this too b/c need modrange(), at least
  %jaaba_modpath_script_path = fullfile(jaaba_perframe_folder_path, 'SetUpJAABAPath.m') ;
  %run(jaaba_modpath_script_path) ;    
  
  % Add FlyBowlAnalysis subfolders that are not their own projects
  addpath(fullfile(path_to_this_folder, 'simplewing')) ;
  addpath(fullfile(path_to_this_folder, 'hmm')) ;  
  
  % Finally, add this folder itself
  addpath(path_to_this_folder) ;  % do this so that we don't have to stay in this folder
  % Also useful b/c this file may be called from elsewhere
  
%   % For debugging, run the JAABA modpath-like script again, b/c DetectJAABA() oes
%   % this anyway
%   jaaba_perframe_folder_path = fullfile(path_to_this_folder, 'JAABA', 'perframe') ;
%   jaaba_modpath_script_path = fullfile(jaaba_perframe_folder_path, 'SetUpJAABAPath.m') ;
%   run(jaaba_modpath_script_path) ;    
end
