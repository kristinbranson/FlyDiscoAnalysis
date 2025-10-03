function success = FlyDiscoAutomaticChecksComplete(expdir, varargin)

version = '0.1';
timestamp = datestr(now,'yyyymmddTHHMMSS');

success = true ;  % Reflects whether any ACC failures have been detected yet (missing *required* files count, missing *desired* files do not)
error_or_warning_messages = {} ;  % Each message in this list is either an error message or a warning message

[analysis_protocol,settingsdir,datalocparamsfilestr,debug_acc,required_file_names_from_stage_name,do_run,~,~] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir',default_settings_folder_path(),...
  'datalocparamsfilestr','dataloc_params.txt',...
  'debug_acc',false, ...    % special debug variable for ACC stage, not the pipeline-wide debug variable
  'required_file_names_from_stage_name', [], ...
  'do_run', [], ...
  'debug', false, ...
  'forcecompute', false);

if ischar(debug_acc),
  debug_acc = (str2double(debug_acc) ~= 0) ;
end
if isempty(required_file_names_from_stage_name) ,
  error('Must supply required_file_names_from_stage_name argument') ;
end
if isempty(do_run) ,
  error('Must supply do_run argument') ;
end

%% parameters
analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
datalocparamsfile = fullfile(analysis_protocol_folder_path, datalocparamsfilestr) ;
dataloc_params = ReadParams(datalocparamsfile) ;
paramsfile = fullfile(analysis_protocol_folder_path, dataloc_params.automaticcheckscompleteparamsfilestr) ;
check_params = ReadParams(paramsfile);
% These fields are not longer used, so delete them
check_params = rmfield(check_params, 'required_files') ;
check_params = rmfield(check_params, 'file_categories') ;
%metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
automatedchecksincomingfile = fullfile(expdir,dataloc_params.automaticchecksincomingresultsfilestr);
acc_text_output_file_path = fullfile(expdir,dataloc_params.automaticcheckscompleteresultsfilestr);

% There's a category for each stage
stage_name_from_stage_index = FlyDiscoStageNames() ;
stage_count = numel(stage_name_from_stage_index) ;
some_categories = cell(1,0) ;
for stage_index = 1 : stage_count ,
  stage_name = stage_name_from_stage_index{stage_index} ;
  some_categories{end+1} = sprintf('missing_%s_files', stage_name) ;  %#ok<AGROW> 
end
% types of automated errors: order matters in this list
other_categories = {...
  'flytracker_nans',...
  'bad_number_of_flies',...
  'bad_number_of_flies_per_sex',...
  'bad_number_of_stimuli',...
  'completed_checks_other'};
categories = horzcat(some_categories, other_categories) ;
category2idx = struct() ;
for i = 1:numel(categories),
  category2idx.(categories{i}) = i;
end
iserror = false(1,numel(categories));

% % make sure there is a category for each file
% if numel(check_params.required_files) ~= numel(check_params.file_categories),
%   error('required_files and file_categories parameters do not match');
% end



%% do some sex-classification-related checks, if that stage is turned on in the screen
if is_on_or_force(do_run.sexclassification) ,
  sexclassifierfile = fullfile(expdir,dataloc_params.sexclassifierdiagnosticsfilestr);
  if ~exist(sexclassifierfile,'file'),
    error_or_warning_messages{end+1} = sprintf('sex classifier diagnostics file %s does not exist',sexclassifierfile);
    success = false;
    iserror(category2idx.missing_sexclassification_files) = true;
  else
    sexclassifier_diagnostics = ReadParams(sexclassifierfile);
    if ~isfield(sexclassifier_diagnostics,'mean_nflies'),
      error_or_warning_messages{end+1} = sprintf('sex classifier diagnostics file %s missing field mean_nflies',sexclassifierfile);
      success = false;
      iserror(category2idx.completed_checks_other) = true;
    else
      num_flies = round(sexclassifier_diagnostics.mean_nflies);
      if num_flies < check_params.min_num_flies,
        error_or_warning_messages{end+1} = sprintf('num_flies = round(%f) = %d < %d',...
          sexclassifier_diagnostics.mean_nflies,num_flies,check_params.min_num_flies);
        success = false;
        iserror(category2idx.bad_number_of_flies) = true;
      end
      if num_flies > check_params.max_num_flies,
        error_or_warning_messages{end+1} = sprintf('num_flies = round(%f) = %d > %d',...
          sexclassifier_diagnostics.mean_nflies,num_flies,check_params.max_num_flies);
        success = false;
        iserror(category2idx.bad_number_of_flies) = true;
      end
    end

    if ~isfield(sexclassifier_diagnostics,'mean_nfemales'),
      error_or_warning_messages{end+1} = sprintf('sex classifier diagnostics file %s missing field mean_nfemales',sexclassifierfile);
      success = false;
      iserror(category2idx.completed_checks_other) = true;
    else
      num_females = round(sexclassifier_diagnostics.mean_nfemales);
      if num_females < check_params.min_num_females,
        error_or_warning_messages{end+1} = sprintf('num_females = round(%f) = %d < %d',...
          sexclassifier_diagnostics.mean_nfemales,num_females,check_params.min_num_females);
        success = false;
        iserror(category2idx.bad_number_of_flies_per_sex) = true;
      end
    end

    if ~isfield(sexclassifier_diagnostics,'mean_nmales'),
      error_or_warning_messages{end+1} = sprintf('sex classifier diagnostics file %s missing field mean_nmales',sexclassifierfile);
      success = false;
      iserror(category2idx.completed_checks_other) = true;
    else
      num_males = round(sexclassifier_diagnostics.mean_nmales);
      if num_males < check_params.min_num_males,
        error_or_warning_messages{end+1} = sprintf('num_males = round(%f) = %d < %d',...
          sexclassifier_diagnostics.mean_nmales,num_males,check_params.min_num_males);
        success = false;
        iserror(category2idx.bad_number_of_flies_per_sex) = true;
      end
    end
  end
end



%% check for nan's in flytracker data after postprocessing 

if is_on_or_force(do_run.computeperframefeatures) ,
  % using the wingtrxfilestr - this is output by perframe features and is the
  % last modification of data in the trx file in the pipeline
  wingtrxfile = fullfile(expdir,dataloc_params.wingtrxfilestr);
  has_wing_info = true ;
  [success, iserror, error_or_warning_messages] = ...
    check_for_nans_in_tracking(wingtrxfile, has_wing_info, category2idx, success, iserror, error_or_warning_messages) ;
elseif is_on_or_force(do_run.registration) ,
  % If PFFs were not calculated in the screen, fallback to checking the registered_trx.mat
  % file.
  registered_trx_file_name = fullfile(expdir,dataloc_params.trxfilestr);
  has_wing_info = false ;
  [success, iserror, error_or_warning_messages] = ...
    check_for_nans_in_tracking(registered_trx_file_name, has_wing_info, category2idx, success, iserror, error_or_warning_messages) ;  
end



%% check the number of stimuli detecter match the expected number from protocol

if is_on_or_force(do_run.ledonoffdetection)
  datalocparamsfile = fullfile(settingsdir,analysis_protocol,'dataloc_params.txt');
  dataloc_params = ReadParams(datalocparamsfile);
  indicatorparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.indicatorparamsfilestr);

  % Get one thing from the indicator params
  if exist(indicatorparamsfile,'file'),
    indicator_params = loadIndicatorParams(indicatorparamsfile) ;
    isOptogeneticExp = logical(indicator_params.OptogeneticExp) ;
  else
    isOptogeneticExp = false ;
    % add error here too?
  end
  if isOptogeneticExp
    [do_stim_counts_match, message] = does_protocol_pulse_count_equal_measured(expdir, settingsdir, analysis_protocol) ;
    if ~do_stim_counts_match ,
      error_or_warning_messages{end+1} = message ;
      success = false;
      iserror(category2idx.bad_number_of_stimuli) = true;
    end
  end
end



%% check for missing files

for i = 1:stage_count ,
  stage_name = stage_name_from_stage_index{i} ; 
  if strcmp(stage_name, 'automaticcheckscomplete') ,
    % These files will typically be missing, b/c we're in the middle
    % of the stage right now!
    continue
  end
  if strcmp(stage_name, 'cleanup') ,
    % This stage runs after ACC.  (And also generates no new files.)
    continue
  end  
  if is_on_or_force(do_run.(stage_name)) ,
    required_file_names = required_file_names_from_stage_name.(stage_name) ;
    for j = 1 : numel(required_file_names) ,
      fn = required_file_names{j} ;
      isfile = logical(exist(fullfile(expdir,fn), 'file')) ;
      if ~isfile ,
        error_or_warning_messages{end+1} = sprintf('Missing file %s',fn) ;  %#ok<AGROW>
        success = false ;
        category = sprintf('missing_%s_files',stage_name) ;
        iserror(category2idx.(category)) = true ;
      end
    end
  end
end



%% check for missing desired files

% To heck with desired files
% if isfield(check_params,'desired_files'),
%   for i = 1:numel(check_params.desired_files),
%     fn = check_params.desired_files{i};
%     if any(fn == '*'),
%       isfile = ~isempty(dir(fullfile(expdir,fn)));
%     else
%       isfile = exist(fullfile(expdir,fn),'file');
%     end
%     if ~isfile,
%       error_or_warning_messages{end+1} = sprintf('Missing desired file %s',fn); %#ok<AGROW>
%       % not that we don't set success = false, b/c these are desired files, not
%       % required files
%     end
%   end
% end



%% output results to file, merging with automatedchecks incoming

if exist(automatedchecksincomingfile,'file'),
  automatedchecks_incoming = ReadParams(automatedchecksincomingfile);
else
  fprintf('Automated checks incoming file %s not find, assuming automated_pf incoming = U\n',automatedchecksincomingfile);
  automatedchecks_incoming = struct('automated_pf','U','notes_curation','');
end

if debug_acc,
  fid = 1;
else
  if exist(acc_text_output_file_path,'file'),
    try %#ok<TRYNC>
      delete(acc_text_output_file_path);
    end
  end
  fid = fopen(acc_text_output_file_path,'w');
end
if fid < 0,
  error('Could not open automatic checks results file %s for writing.',acc_text_output_file_path);
end
if success && ~strcmpi(automatedchecks_incoming.automated_pf,'F'),
  fprintf(fid,'automated_pf,P\n');
  if ~isempty(error_or_warning_messages)
      fprintf(fid,'notes_curation,');
      s = sprintf('%s\\n',error_or_warning_messages{:});
      s = s(1:end-2);
      if isfield(automatedchecks_incoming,'notes_curation') && ...
              ~isempty(automatedchecks_incoming.notes_curation),
          if isempty(s),
              s = automatedchecks_incoming.notes_curation;
          else
              s = [automatedchecks_incoming.notes_curation,'\n',s];
          end
      end
      fprintf(fid,'%s\n',s);
  end
else
  fprintf(fid,'automated_pf,F\n');
  fprintf(fid,'notes_curation,\n');
  s = sprintf('%s,\\n',error_or_warning_messages{:});
%   s = s(1:end-2);
  if isfield(automatedchecks_incoming,'notes_curation') && ...
      ~isempty(automatedchecks_incoming.notes_curation),
    if isempty(s),
      s = automatedchecks_incoming.notes_curation;
    else
      s = [automatedchecks_incoming.notes_curation,'\n',s];
    end
  end
% AR 20240530 changed autochecks complete formatting to make it easier to read 
  fprintf(fid,s);
  if strcmpi(automatedchecks_incoming.automated_pf,'F'),
    if isfield(automatedchecks_incoming,'automated_pf_category') && ...
        ~isempty(automatedchecks_incoming.automated_pf_category),
      s = automatedchecks_incoming.automated_pf_category;
    else
      s = 'incoming_checks_unknown';
    end
  else
    i = find(iserror,1);
    if isempty(i),
      s = 'completed_checks_other';
    else
      s = categories{i};
    end
  end
  fprintf(fid,'automated_pf_category,%s\n',s);      
  if debug_acc,
    fprintf('[[automated_pf_categories:');
    fprintf(' %s',categories{iserror});
    fprintf(']]\n');
  end
end

if ~debug_acc,
  fclose(fid);
end



%% save info to mat file

real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);
acc_mat_output_file_path = fullfile(expdir,dataloc_params.automaticcheckscompleteinfomatfilestr);
fprintf('Saving debug info to file %s...\n',acc_mat_output_file_path);

accinfo = struct() ;
accinfo.paramsfile = paramsfile;
accinfo.check_params = check_params;
accinfo.version = version;
accinfo.analysis_protocol = analysis_protocol;
accinfo.linked_analysis_protocol = real_analysis_protocol;
accinfo.timestamp = timestamp;
accinfo.iserror = iserror;
accinfo.categories = categories;
accinfo.msgs = error_or_warning_messages;
accinfo.automatedchecks_incoming = automatedchecks_incoming;
accinfo.success = success ;

if ~debug_acc,
  if exist(acc_mat_output_file_path,'file'),
    try %#ok<TRYNC>
      delete(acc_mat_output_file_path);
    end
  end
  try
    save(acc_mat_output_file_path,'-struct','accinfo');
  catch ME,
    warning('Could not save information to file %s: %s',acc_mat_output_file_path,getReport(ME));
  end
end



%% print results to log file

%fprintf('success = %d\n',success);
stage = 'automaticchecks_complete' ;
fprintf('Stage %s success = %d\n', stage, success) ;
if success ,
  if isempty(error_or_warning_messages) ,     
    fprintf('%s: No error or warning messages.\n', stage);
  else
    fprintf('%s: No errors, but there were warnings:\n', stage);
    fprintf('%s\n',error_or_warning_messages{:});      
  end
else
  %flydisco_pipeline_error(stage, error_or_warning_messages) ;
  fprintf('Warning/error messages:\n');
  fprintf('%s\n',error_or_warning_messages{:});
end
fprintf('Finished running FlyDiscoAutomaticChecksComplete at %s.\n',datestr(now,'yyyymmddTHHMMSS'));

