function modpath()
  % Add needed libraries to Matlab path
  
  path_to_this_script = mfilename('fullpath') ;
  path_to_this_folder = fileparts(path_to_this_script) ;
  
%   % Run the FlyBowlAnalysis modpath script
%   fly_bowl_analysis_folder_path = fullfile(path_to_this_folder, 'FlyDiscoAnalysis') ;
%   fly_bowl_analysis_modpath_script_path = fullfile(fly_bowl_analysis_folder_path, 'modpath.m') ;
%   run(fly_bowl_analysis_modpath_script_path) ;  
  
  % Add the fuster folder
  addpath(fullfile(path_to_this_folder, 'fuster')) ;
  
  % Add this folder
  addpath(path_to_this_folder) ;  % do this so that we don't have to stay in this folder
  % Also useful b/c this file may be called from elsewhere
  