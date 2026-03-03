function FlyDiscoCleanup(expdir, varargin)
% Stage to cleanup files we don't want to keep long-term.  Runs only if the ACC
% stage reports the pipeline ran cleanly.

% Parse optional args
[analysis_protocol, settingsdir, datalocparamsfilestr, ~, ~, do_run, ~, ~] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir',default_settings_folder_path(),...
  'datalocparamsfilestr','dataloc_params.txt',...
  'debug_acc',false, ...    % special debug variable for ACC stage, not the pipeline-wide debug variable
  'required_file_names_from_stage_name', [], ...
  'do_run', [], ...
  'debug', false, ...
  'forcecompute', false);

% Read the stage parameters
analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
datalocparamsfile = fullfile(analysis_protocol_folder_path, datalocparamsfilestr) ;
dataloc_params = ReadParams(datalocparamsfile) ;

% If dataloc_params doesn't specify the cleanup params file, return
if ~isfield(dataloc_params, 'cleanupparamsfilestr')
  % Seems cleanup stage params are not defined in dataloc_params, so presumably
  % this analysis protocol predates the cleanup stage, or users don't need any
  % cleanup.
  fprintf('The file %s does not define cleanupparamsfilestr, so not doing any cleanup.\n', datalocparamsfilestr) ;
  return
end

% Check if we should proceed with cleanup.
% Checks if the ACC stage is on, the results file is present and indicates the pipeline passed.
do_run_cleanup_stage = CheckACCFile(expdir, dataloc_params, do_run) ;
if ~do_run_cleanup_stage
  fprintf('The ACC stage is off, or the ACC result file does not exist, or does not indicate the pipeline ran successfully, so not doing any cleanup.\n') ;
  return  
end

% Read the stage parameters
cleanup_params_file_name = dataloc_params.cleanupparamsfilestr ;
cleanup_params_file_path = fullfile(analysis_protocol_folder_path, cleanup_params_file_name) ;
if ~exist(cleanup_params_file_path, 'file')
  fprintf('The cleanup stage parameter file %s does not exist, so not doing any cleanup.\n', cleanup_params_file_name) ;
  return
end
stage_params = ReadParams(cleanup_params_file_path);

% Generate the black list
if isfield(stage_params, 'blacklist_glob_file_name')
  blacklist_glob_file_name = stage_params.blacklist_glob_file_name ;  % the name of the file containing the list of globs
  blacklist_file_paths = path_list_from_color_file_name(blacklist_glob_file_name, analysis_protocol_folder_path, expdir) ;
else
  blacklist_file_paths = cell(1,0) ;
end

% Generate the white list
if isfield(stage_params, 'whitelist_glob_file_name')
  whitelist_glob_file_name = stage_params.whitelist_glob_file_name ;  % the name of the file containing the list of globs
  whitelist_file_paths = path_list_from_color_file_name(whitelist_glob_file_name, analysis_protocol_folder_path, expdir) ;
else
  whitelist_file_paths = cell(1,0) ;
end

% The white list takes priority, so remove any whitelist items in the
% blacklist
delete_list_file_paths = setdiff(blacklist_file_paths, whitelist_file_paths) ;

% Delete the files
to_delete_count = numel(delete_list_file_paths) ;
deleted_count = 0 ;
for  i = 1 : to_delete_count
  file_path = delete_list_file_paths{i} ;
  try
    delete(file_path) ;
    deleted_count = deleted_count + 1 ;
  catch me
    fprintf('During cleanup, unable to delete file %s: %s', file_path, me.msg) ;
  end
end
if deleted_count > 0 
  fprintf('Deleted %d files.\n', deleted_count) ;
end

end  % function
