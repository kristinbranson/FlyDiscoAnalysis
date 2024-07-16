function FlyDiscoJAABADetect(expdir, varargin)

% Process optional args
[analysis_protocol, settingsdir, forcecompute] = ...
  myparse(varargin,...
          'analysis_protocol', '', ...
          'settingsdir', default_settings_folder_path(), ...
          'forcecompute', false) ;
if isempty(analysis_protocol) ,
  error('analysis_protocol cannot be empty') ;
end

% Read in the names of the .jab files
datalocparamsfilestr = 'dataloc_params.txt' ;
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr) ;
dataloc_params = ReadParams(datalocparamsfile) ;
jaabaclassifierparamsfilestrs = fullfile(settingsdir,analysis_protocol,dataloc_params.jaabaclassifierparamsfilestrs) ;
jabfiles = read_one_file_name_per_line(jaabaclassifierparamsfilestrs) ;
  
% Actually call JAABADetect()
% For reasons that are unclear to me, the otherwise-useless try/catch
% wrapper ensures that errors encounted during loading of user-defined
% class objects are still silently ignored, even when "dbstop if error"
% is engaged.  ALT, 2021-03-04
try
  JAABADetect(expdir, 'jabfiles', jabfiles, 'forcecompute', forcecompute) ;
catch me ,
  rethrow(me) ;
end

end  % function
