function FlyDiscoJAABADetect(expdir, varargin)

% Process optional args
[analysis_protocol, settingsdir, forcecompute, ~, ~] = ...
  myparse(varargin,...
          'analysis_protocol', '', ...
          'settingsdir', default_settings_folder_path(), ...
          'forcecompute', false, ...
          'debug', false, ...
          'do_run', []) ;
if isempty(analysis_protocol) ,
  error('analysis_protocol cannot be empty') ;
end

% Read the data location parameters
datalocparamsfilestr = 'dataloc_params.txt' ;
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr) ;
dataloc_params = ReadParams(datalocparamsfile) ;

% Read in the names of the .jab files
jaabaclassifierparamsfilestrs = fullfile(settingsdir,analysis_protocol,dataloc_params.jaabaclassifierparamsfilestrs) ;
jabfiles = read_one_file_name_per_line(jaabaclassifierparamsfilestrs) ;

% Set the fastcomputepffs option
if isfield(dataloc_params, 'jaabadetectparamsfilestr')
  paramsfile = fullfile(settingsdir, analysis_protocol, dataloc_params.jaabadetectparamsfilestr) ;
  if exist(paramsfile, 'file')
    params = ReadParams(paramsfile) ;
    if isfield(params, 'fastcomputepffs')
      fastcomputepffs = logical(params.fastcomputepffs) ;
    else
      fastcomputepffs = false ;
    end
  else
    fastcomputepffs = false ;
  end
else
  fastcomputepffs = false ;
end

% Actually call JAABADetect()
% For reasons that are unclear to me, the otherwise-useless try/catch
% wrapper ensures that errors encounted during loading of user-defined
% class objects are still silently ignored, even when "dbstop if error"
% is engaged.  ALT, 2021-03-04
try
  JAABADetect(expdir, 'jabfiles', jabfiles, 'forcecompute', forcecompute, 'fastcomputepffs', fastcomputepffs) ;
catch me ,
  rethrow(me) ;
end

end  % function
