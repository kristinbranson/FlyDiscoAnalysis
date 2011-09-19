function [success,msgs] = FlyBowlAutomaticChecks_Complete(expdir,varargin)

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
paramsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.automaticcheckscompleteparamsfilestr);
check_params = ReadParams(paramsfile);
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
automatedchecksincomingfile = fullfile(expdir,dataloc_params.automaticchecksincomingresultsfilestr);
outfile = fullfile(expdir,dataloc_params.automaticcheckscompleteresultsfilestr);

% types of automated errors: order matters in this list
categories = {...
  'missing_automated_checks_incoming_files',...
  'missing_tracking_files',...
  'ctrax_infinity_bug',...
  'missing_registration_files',...
  'missing_sexclassification_files',...
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

if ~isfield(metadata,'screen_type'),
  success = false;
  msgs{end+1} = 'screen_type not stored in Metadata file';
  isscreen = false;
else
  isscreen = ~strcmpi(metadata.screen_type,'non_olympiad') && ...
    ~strcmpi(metadata.screen_type,'non_production');
end

%% check number of flies

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

%% check for Ctrax bugs: nan, inf in trajectories

ctraxfile = fullfile(expdir,dataloc_params.ctraxfilestr);
if ~exist(ctraxfile,'file'),
  msgs{end+1} = sprintf('Ctrax output mat file %s does not exist',ctraxfile);
  success = false;
  iserror(category2idx.missing_tracking_files) = true;
else
  % name of annotation file
  annfile = fullfile(expdir,dataloc_params.annfilestr);
  
  % name of movie file
  moviefile = fullfile(expdir,dataloc_params.moviefilestr);
  
  % load trajectories
  [trx,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile,'annname',annfile); %#ok<NASGU>
  if ~succeeded,
    msgs{end+1} = sprintf('Could not load trajectories from file %s',ctraxfile);
    success = false;
    iserror(category2idx.completed_checks_other) = true;
  else
    
    for fly = 1:numel(trx),
      badidx = isnan(trx(fly).x) | ...
        isnan(trx(fly).y) | ...
        isnan(trx(fly).a) | ...
        isnan(trx(fly).b) | ...
        isnan(trx(fly).theta);
      if any(badidx),
        [starts,ends] = get_interval_ends(badidx);
        starts = starts - trx(fly).off;
        ends = ends - trx(fly).off;
        msgs{end+1} = [sprintf('Trajectory %d has NaNs in frames',fly),sprintf(' %d-%d',[starts,ends]')]; %#ok<AGROW>
        success = false;
        iserror(category2idx.ctrax_infinity_bug) = true;
      end
      badidx = isinf(trx(fly).x) | ...
        isinf(trx(fly).y) | ...
        isinf(trx(fly).a) | ...
        isinf(trx(fly).b) | ...
        isinf(trx(fly).theta);
      if any(badidx),
        [starts,ends] = get_interval_ends(badidx);
        starts = starts - trx(fly).off;
        ends = ends - trx(fly).off;
        msgs{end+1} = [sprintf('Trajectory %d has Infs in frames',fly),sprintf(' %d-%d',[starts,ends]')]; %#ok<AGROW>
        success = false;
        iserror(category2idx.ctrax_infinity_bug) = true;
      end
      
    end
  end
end

%% check for missing files

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

%% output results to file, merging with automatedchecks incoming

if exist(automatedchecksincomingfile,'file'),
  automatedchecks_incoming = ReadParams(automatedchecksincomingfile);
else
  automatedchecks_incoming = struct('automated_pf','P','notes_curation','');
end

if DEBUG,
  fid = 1;
else
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
  if strcmpi(automatedchecks_incoming.automated_pf,'F') && ...
      isfield(automatedchecks_incoming,'automated_pf_category') && ...
      ~isempty(automatedchecks_incoming.automated_pf_category),
    s = automatedchecks_incoming.automated_pf_category;
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
    fprintf('automated_pf_categories:');
    fprintf(' %s',categories{iserror});
    fprintf('\n');
  end
end

if ~DEBUG,
  fclose(fid);
end