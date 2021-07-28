function [success,msgs,iserror] = FlyDiscoAutomaticChecksComplete(expdir,varargin)

version = '0.1';
timestamp = datestr(now,'yyyymmddTHHMMSS');



success = true;
msgs = {};

[analysis_protocol,settingsdir,datalocparamsfilestr,DEBUG] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt','debug',false);

if ischar(DEBUG),
  DEBUG = str2double(DEBUG) ~= 0;
end

%% parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% log file
% 
% if isfield(dataloc_params,'automaticchecks_complete_logfilestr') && ~DEBUG,
%   logfile = fullfile(expdir,dataloc_params.automaticchecks_complete_logfilestr);
%   logfid = fopen(logfile,'a');
%   if logfid < 1,
%     warning('Could not open log file %s\n',logfile);
%     logfid = 1;
%   end
% else
%   logfid = 1;
% end
% 
real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);
% 
% fprintf(logfid,'\n\n***\nRunning FlyDiscoAutomaticChecks_Complete version %s analysis_protocol %s (linked to %s) at %s\n',version,analysis_protocol,real_analysis_protocol,timestamp);

  
%% more parameters

paramsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.automaticcheckscompleteparamsfilestr);
check_params = ReadParams(paramsfile);
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
automatedchecksincomingfile = fullfile(expdir,dataloc_params.automaticchecksincomingresultsfilestr);
outfile = fullfile(expdir,dataloc_params.automaticcheckscompleteresultsfilestr);

% types of automated errors: order matters in this list
categories = {...
  'missing_automated_checks_incoming_files',...
  'missing_tracking_files',...
  'flytracker_nans',...
  'missing_registration_files',...
  'missing_optoregistration_files',...
  'missing_sexclassification_files',...
  'missing_wingtracking_files',...
  'missing_perframefeatures_files',...
  'missing_perframestats_files',...
  'missing_extra_diagnostics_files',...
  'missing_results_movie_files',...
  'missing_other_analysis_files',...
  'bad_number_of_flies',...
  'bad_number_of_flies_per_sex',...
  'completed_checks_other'};
category2idx = struct;
for i = 1:numel(categories),
  category2idx.(categories{i}) = i;
end
iserror = false(1,numel(categories));

% make sure there is a category for each file
if numel(check_params.required_files) ~= numel(check_params.file_categories),
  error('required_files and file_categories parameters do not match');
end

%% read metadata

metadata = ReadMetadataFile(metadatafile);

%% check for primary screen: some checks are only for screen data

% if ~isfield(metadata,'screen_type'),
%   success = false;
%   msgs{end+1} = 'screen_type not stored in Metadata file';
%   isscreen = false;
% else
%   isscreen = ~strcmpi(metadata.screen_type,'non_olympiad') && ...
%     ~strcmpi(metadata.screen_type,'non_production');
% end
% TODO change this as noted below 
isscreen = true;
%% check number of flies if flies nums included in autochecks 
% TODO edit this to check if numflies values are included in
% autocheckscomplete params. 

if isscreen,
  
  sexclassifierfile = fullfile(expdir,dataloc_params.sexclassifierdiagnosticsfilestr);
  if ~exist(sexclassifierfile,'file'),
    msgs{end+1} = sprintf('sex classifier diagnostics file %s does not exist',sexclassifierfile);
    success = false;
    iserror(category2idx.missing_sexclassification_files) = true;
  else
    sexclassifier_diagnostics = ReadParams(sexclassifierfile);
    if ~isfield(sexclassifier_diagnostics,'mean_nflies'),
      msgs{end+1} = sprintf('sex classifier diagnostics file %s missing field mean_nflies',sexclassifierfile);
      success = false;
      iserror(category2idx.completed_checks_other) = true;
    else
      num_flies = round(sexclassifier_diagnostics.mean_nflies);
      if num_flies < check_params.min_num_flies,
        msgs{end+1} = sprintf('num_flies = round(%f) = %d < %d',...
          sexclassifier_diagnostics.mean_nflies,num_flies,check_params.min_num_flies);
        success = false;
        iserror(category2idx.bad_number_of_flies) = true;
      end
      if num_flies > check_params.max_num_flies,
        msgs{end+1} = sprintf('num_flies = round(%f) = %d > %d',...
          sexclassifier_diagnostics.mean_nflies,num_flies,check_params.max_num_flies);
        success = false;
        iserror(category2idx.bad_number_of_flies) = true;
      end
    end
    
    if ~isfield(sexclassifier_diagnostics,'mean_nfemales'),
      msgs{end+1} = sprintf('sex classifier diagnostics file %s missing field mean_nfemales',sexclassifierfile);
      success = false;
      iserror(category2idx.completed_checks_other) = true;
    else
      num_females = round(sexclassifier_diagnostics.mean_nfemales);
      if num_females < check_params.min_num_females,
        msgs{end+1} = sprintf('num_females = round(%f) = %d < %d',...
          sexclassifier_diagnostics.mean_nfemales,num_females,check_params.min_num_females);
        success = false;
        iserror(category2idx.bad_number_of_flies_per_sex) = true;
      end
    end
    
    if ~isfield(sexclassifier_diagnostics,'mean_nmales'),
      msgs{end+1} = sprintf('sex classifier diagnostics file %s missing field mean_nmales',sexclassifierfile);
      success = false;
      iserror(category2idx.completed_checks_other) = true;
    else
      num_males = round(sexclassifier_diagnostics.mean_nmales);
      if num_males < check_params.min_num_males,
        msgs{end+1} = sprintf('num_males = round(%f) = %d < %d',...
          sexclassifier_diagnostics.mean_nmales,num_males,check_params.min_num_males);
        success = false;
        iserror(category2idx.bad_number_of_flies_per_sex) = true;
      end
    end
    
  end
  
end

%% check for nan in flytracker data after postprocessing 

% using the wingtrxfilestr - this is output by perframe features and is the
% last modification of data in the trx file in the pipeline
wingtrxfile = fullfile(expdir,dataloc_params.wingtrxfilestr);
if ~exist(wingtrxfile,'file'),
  msgs{end+1} = sprintf('wingtracking_results mat file %s does not exist',wingtrxfile);
  success = false;
  iserror(category2idx.missing_tracking_files) = true;
else
  % name of annotation file
  %annfile = fullfile(expdir,dataloc_params.annfilestr);
  
  % name of movie file
  %moviefile = fullfile(expdir,dataloc_params.moviefilestr);
  
  % load trajectories
  [trx,~,succeeded,timestamps] = load_tracks(wingtrxfile); %#ok<NASGU>
  if ~succeeded,
    msgs{end+1} = sprintf('Could not load trajectories from file %s',wingtrxfile);
    success = false;
    iserror(category2idx.completed_checks_other) = true;
  else
    
    for fly = 1:numel(trx),
      badidx = isnan(trx(fly).x) | ...
        isnan(trx(fly).y) | ...
        isnan(trx(fly).a) | ...
        isnan(trx(fly).b) | ...
        isnan(trx(fly).theta) | ...
        isnan(trx(fly).xwingl) | ...
        isnan(trx(fly).ywingl) | ...
        isnan(trx(fly).xwingr) | ...
        isnan(trx(fly).ywingr) | ...
        isnan(trx(fly).wing_anglel) | ...
        isnan(trx(fly).wing_angler);    
    
      if any(badidx),
        [starts,ends] = get_interval_ends(badidx);
        starts = starts - trx(fly).off;
        ends = ends - trx(fly).off;
        for se = 1:numel(starts)
            msgs{end+1} = [sprintf('Trajectory %d has NaNs in frames',fly),sprintf(' %d-%d',[starts(se),ends(se)]')]; %#ok<AGROW>
        end
        success = false;
        iserror(category2idx.flytracker_nans) = true;
      end
%       badidx = isinf(trx(fly).x) | ...
%         isinf(trx(fly).y) | ...
%         isinf(trx(fly).a) | ...
%         isinf(trx(fly).b) | ...
%         isinf(trx(fly).theta);
%       if any(badidx),
%         [starts,ends] = get_interval_ends(badidx);
%         starts = starts - trx(fly).off;
%         ends = ends - trx(fly).off;
%         msgs{end+1} = [sprintf('Trajectory %d has Infs in frames',fly),sprintf(' %d-%d',[starts,ends]')]; %#ok<AGROW>
%         success = false;
%         iserror(category2idx.ctrax_infinity_bug) = true;
%       end
      
    end
  end
end

%% check for missing files
% add parsing for opto or not

for i = 1:numel(check_params.required_files),
  fn = check_params.required_files{i};
  if any(fn == '*'),
    isfile = ~isempty(dir(fullfile(expdir,fn)));
  else
    isfile = exist(fullfile(expdir,fn),'file');
  end
  if ~isfile,
    msgs{end+1} = sprintf('Missing file %s',fn); %#ok<AGROW>
    success = false;
    category = sprintf('missing_%s_files',check_params.file_categories{i});
    iserror(category2idx.(category)) = true;
  end
end

%% check for scores files
% 
% try
%   if isfield(dataloc_params,'jaabadetectparamsfilestr'),
%     jaabadetectparams = ReadParams(fullfile(settingsdir,analysis_protocol,dataloc_params.jaabadetectparamsfilestr));
%     [issuccess_scores,msgs_scores] = CheckScores({expdir},fullfile(settingsdir,analysis_protocol,jaabadetectparams.classifierparamsfiles));
%     if any(~issuccess_scores),
%       msgcurr = sprintf('Bad JAABA scores files\n');
%       for i = 1:numel(msgs_scores),
%         if ~isempty(msgs_scores{i}),
%           msgcurr = [msgcurr,sprintf('%s\n',msgs_scores{i}{:})]; %#ok<AGROW>
%         end
%       end
%       msgs{end+1} = msgcurr;
%       success = false;
%       iserror(category2idx.bad_jaaba_scores) = true;
%     elseif ~isempty(msgs_scores),
%       msgs(end+1:end+numel(msgs_scores)) = msgs_scores;
%     end
%   end
% catch ME,
%   msgs{end+1} = sprintf('Error checking JAABA scores files: %s',getReport(ME));
%   success = false;
%   iserror(category2idx.completed_checks_other) = true;
% end

%% check for missing desired files

if isfield(check_params,'desired_files'),
  for i = 1:numel(check_params.desired_files),
    fn = check_params.desired_files{i};
    if any(fn == '*'),
      isfile = ~isempty(dir(fullfile(expdir,fn)));
    else
      isfile = exist(fullfile(expdir,fn),'file');
    end
    if ~isfile,
      msgs{end+1} = sprintf('Missing desired file %s',fn); %#ok<AGROW>
    end
  end
end

%% output results to file, merging with automatedchecks incoming

if exist(automatedchecksincomingfile,'file'),
  automatedchecks_incoming = ReadParams(automatedchecksincomingfile);
else
  fprintf('Automated checks incoming file %s not find, assuming automated_pf incoming = U\n',automatedchecksincomingfile);
  automatedchecks_incoming = struct('automated_pf','U','notes_curation','');
end

if DEBUG,
  fid = 1;
else
  if exist(outfile,'file'),
    try %#ok<TRYNC>
      delete(outfile);
    end
  end
  fid = fopen(outfile,'w');
end
if fid < 0,
  error('Could not open automatic checks results file %s for writing.',outfile);
end
if success && ~strcmpi(automatedchecks_incoming.automated_pf,'F'),
  fprintf(fid,'automated_pf,P\n');
else
  fprintf(fid,'automated_pf,F\n');
  fprintf(fid,'notes_curation,');
  s = sprintf('%s\\n',msgs{:});
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
  if DEBUG,
    fprintf('[[automated_pf_categories:');
    fprintf(' %s',categories{iserror});
    fprintf(']]\n');
  end
end

if ~DEBUG,
  fclose(fid);
end

%% save info to mat file

filename = fullfile(expdir,dataloc_params.automaticcheckscompleteinfomatfilestr);
fprintf('Saving debug info to file %s...\n',filename);

accinfo = struct;
accinfo.paramsfile = paramsfile;
accinfo.check_params = check_params;
accinfo.version = version;
accinfo.analysis_protocol = analysis_protocol;
accinfo.linked_analysis_protocol = real_analysis_protocol;
accinfo.timestamp = timestamp;
accinfo.iserror = iserror;
accinfo.categories = categories;
accinfo.msgs = msgs;
accinfo.automatedchecks_incoming = automatedchecks_incoming;
accinfo.success = success; %#ok<STRNU>

if ~DEBUG,
  if exist(filename,'file'),
    try %#ok<TRYNC>
      delete(filename);
    end
  end
  try
    save(filename,'-struct','accinfo');
  catch ME,
    warning('Could not save information to file %s: %s',filename,getReport(ME));
  end
end

%% print results to log file

fprintf('success = %d\n',success);
if isempty(msgs),
  fprintf('No error or warning messages.\n');
else
  fprintf('Warning/error messages:\n');
  fprintf('%s\n',msgs{:});
end
fprintf('Finished running FlyDiscoAutomaticChecks_Complete at %s.\n',datestr(now,'yyyymmddTHHMMSS'));

