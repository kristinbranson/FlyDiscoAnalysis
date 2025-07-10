function modpath()
  % Add needed libraries to Matlab path

  % Sort out where FlyDiscoAnalysis is in the filesytem
  path_to_this_script = mfilename('fullpath') ;
  path_to_this_folder = fileparts(path_to_this_script) ;
  
%   % Add the fuster folder
%   addpath(fullfile(path_to_this_folder, 'fuster')) ;
  
  % Add TMT to the path, so we get fuster stuff, but do it early so it gets
  % shadowed by anything with a namespace conflict.
  tmt_folder_path = fullfile(path_to_this_folder, 'tmt') ;
  tmt_modpath_script_path = fullfile(tmt_folder_path, 'modpath.m') ; 
  run(tmt_modpath_script_path) ;  

  % Run the FlyTracker modpath script
  flytracker_folder_path = fullfile(path_to_this_folder, 'FlyTracker') ;
  flytracker_modpath_script_path = fullfile(flytracker_folder_path, 'modpath.m') ;
  run(flytracker_modpath_script_path) ;  

  % Add the JAABA stuff to the path
  jaaba_perframe_folder_path = fullfile(path_to_this_folder, 'JAABA', 'perframe') ;
  jaaba_modpath_script_path = fullfile(jaaba_perframe_folder_path, 'SetUpJAABAPath.m') ;
  run(jaaba_modpath_script_path) ;    

  % Remove spaceTime stuff...
  w = warning('query','MATLAB:rmpath:DirNotFound');
  warning('off','MATLAB:rmpath:DirNotFound');
  rmpath(genpath(fullfile(path_to_this_folder, 'JAABA', 'spaceTime'))) ;
  warning(w.state,w.identifier);

  % Add the TrkFile code for loading in trk files
  % The -begin option forces it to the front of the path, even if it's on the
  % path already.
  addpath(fullfile(path_to_this_folder,'APT','matlab','trk'), '-begin') ;

  % Add FlyDiscoAnalysis subfolders that are not their own projects
  addpath(fullfile(path_to_this_folder, 'simplewing'), '-begin') ;
  addpath(fullfile(path_to_this_folder, 'hmm'), '-begin') ;  
  addpath(fullfile(path_to_this_folder, 'filehandling'), '-begin') ;  
  addpath(fullfile(path_to_this_folder, 'perframe'), '-begin') ;
  addpath(fullfile(path_to_this_folder, 'utility'), '-begin') ;  
  addpath(fullfile(path_to_this_folder, 'locomotion'), '-begin') ; 
  addpath(fullfile(path_to_this_folder, 'locomotion/external'), '-begin') ; 

  % Add stuff intended to shadow JAABA, FlyTracker versions
  addpath(fullfile(path_to_this_folder, 'shadow'), '-begin') ;
  
  % Add various folders
  addpath(fullfile(path_to_this_folder, 'tests'), '-begin') ;
  addpath(fullfile(path_to_this_folder, 'scripts'), '-begin') ;
  addpath(fullfile(path_to_this_folder, 'debug'), '-begin') ;
  addpath(fullfile(path_to_this_folder, 'main'), '-begin') ;
  addpath(fullfile(path_to_this_folder, 'stages'), '-begin') ;
  
  % Finally, add this folder itself, so we don't have to stay in this folder
  % Add at the beginning so that e.g. FlyTracker doesn't override
  addpath(path_to_this_folder, '-begin') ;
  
  % Run code to set the parpool location appropriately
  set_parpool_job_storage_location()
end
