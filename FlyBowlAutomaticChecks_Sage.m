function [success,msgs] = FlyBowlAutomaticChecks_Sage(expdir,varargin)

success = true;
msgs = {};

[analysis_protocol,settingsdir,datalocparamsfilestr,first_barcode_datetime,dosave] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'first_barcode_datetime','20110413T000000',...
  'dosave',true);

if ischar(dosave),
  dosave = str2double(dosave) ~= 0;
end
datetime_format = 'yyyymmddTHHMMSS';

%% parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
paramsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.automaticchecksincomingparamsfilestr);
check_params = ReadParams(paramsfile);
if ~iscell(check_params.control_line_names),
  check_params.control_line_names = {check_params.control_line_names};
end
%metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
%ufmfdiagnosticsfile = fullfile(expdir,dataloc_params.ufmfdiagnosticsfilestr);
%temperaturefile = fullfile(expdir,dataloc_params.temperaturefilestr);
if dosave,
  outfile_incoming = fullfile(expdir,dataloc_params.automaticchecksincomingresultsfilestr);
  outfile_complete = fullfile(expdir,dataloc_params.automaticcheckscompleteresultsfilestr);
else
  outfile_incoming = '';
  outfile_complete = '';
end

%% read metadata

[~,experiment_name] = fileparts(expdir);
experiment_name = ['FlyBowl_',experiment_name];
metadata = SAGEGetBowlData('dataset','score','experiment_name',experiment_name,'checkflags',false);

%% get notes_curations
if ~isempty(metadata.notes_curation),
  msgs{end+1} = metadata.notes_curation;
end

%% check for metadata fields

required_fns = {'flag_aborted','flag_redo','seconds_fliesloaded','seconds_shiftflytemp',...
  'screen_type','line_name','cross_barcode'};
ismissingfn = ~ismember(required_fns,fieldnames(metadata));
if any(ismissingfn),
  success = false;
  msgs{end+1} = ['Missing metadata fields:',sprintf(' %s',required_fns{ismissingfn})];
end

%% check for flags

if metadata.flag_aborted ~= 0,
  success = false;
  msgs{end+1} = 'Experiment aborted.';
end

if metadata.flag_redo ~= 0,
  success = false;
  msgs{end+1} = 'Redo flag set to 1.';
end

%% check loading time

if metadata.seconds_fliesloaded < check_params.min_seconds_fliesloaded,
  success = false;
  msgs{end+1} = sprintf('Load time = %f < %f seconds.',metadata.seconds_fliesloaded,check_params.min_seconds_fliesloaded);
end
if metadata.seconds_fliesloaded > check_params.max_seconds_fliesloaded,
  success = false;
  msgs{end+1} = sprintf('Load time = %f > %f seconds.',metadata.seconds_fliesloaded,check_params.max_seconds_fliesloaded);
end

%% check shiftflytemp time
if metadata.seconds_shiftflytemp > check_params.max_seconds_shiftflytemp,
  success = false;
  msgs{end+1} = sprintf('Shift fly temp time = %f > %f seconds.',metadata.seconds_shiftflytemp,check_params.max_seconds_shiftflytemp);
end

%% check for primary screen: all other checks are only for screen data

if ~isfield(metadata,'screen_type'),
  success = false;
  msgs{end+1} = 'screen_type not stored in Metadata file';
  isscreen = false;
else
  isscreen = ~strcmpi(metadata.screen_type,'non_olympiad') && ...
    ~strcmpi(metadata.screen_type,'non_production');
end

%% check barcode

if isscreen,

if ~ismember(metadata.line_name,check_params.control_line_names) && ...
    datenum(metadata.exp_datetime,datetime_format) > datenum(first_barcode_datetime,datetime_format) && ...
    str2double(metadata.cross_barcode) < 0,
  success = false;
  msgs{end+1} = 'Barcode = -1 and line_name indicates not control';
end

%% check video length

if ~isfield(metadata,'ufmf_diagnostics_summary_nFrames'),
  success = false;
  msgs{end+1} = 'ufmf_diagnostics_summary_nFrames missing.';
else
  if metadata.ufmf_diagnostics_summary_nFrames < check_params.min_ufmf_diagnostics_summary_nframes,
    success = false;
    msgs{end+1} = sprintf('Video contains %d < %d frames.',metadata.ufmf_diagnostics_summary_nFrames,check_params.min_ufmf_diagnostics_summary_nframes);
  end
end

%% check temperature
% 
% if ~exist(temperaturefile,'file'),
%   warning('Temperature file %s does not exist',temperaturefile);
% else
%   try
%     tempdata = importdata(temperaturefile,',');
%   catch ME,
%     warning('Error importing temperature stream: %s',getReport(ME));
%   end
%   if isempty(tempdata),
%     warning('No temperature readings recorded.');
%     temp = [];
%   elseif size(tempdata,2) < 2,
%     warning('Temperature data could not be read');
%     temp = [];
%   else
%     temp = tempdata(:,2);
%     if isempty(temp),
%       warning('No temperature readings recorded.');
%     else
%       if max(temp) > check_params.max_temp,
%         success = false;
%         msgs{end+1} = sprintf('Max temperature = %f > %f.',max(temp),check_params.max_temp);
%       end
%       if numel(temp) < 2,
%         warning('Only one temperature recorded.');
%       else
%         if max(temp) - min(temp) > check_params.max_tempdiff,
%           success = false;
%           msgs{end+1} = sprintf('Temperature change = %f > %f.',max(temp) - min(temp),check_params.max_tempdiff);
%         end
%       end
%     end
%   end
% end

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
  end
end

%% output results to file

if isempty(outfile_incoming),
  fid = 1;
else
  fid = fopen(outfile_incoming,'w');
end
if fid < 0,
  error('Could not open automatic checks results file %s for writing.',outfile);
end
if success,
  fprintf(fid,'automated_pf,P\n');
else
  fprintf(fid,'automated_pf,F\n');
  fprintf(fid,'notes_curation,');
  s = sprintf('%s\\n',msgs{:});
  s = s(1:end-2);
  fprintf(fid,'%s\n',s);
end

if fid > 1,
  fclose(fid);
end

%% complete checks

paramsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.automaticcheckscompleteparamsfilestr);
check_params = ReadParams(paramsfile);

success_incoming = success;
msgs_incoming = msgs;
success = true;
msgs = {};

%% check number of flies

if isscreen,
  
  if ~isfield(metadata,'sexclassifier_diagnostics_mean_nflies'),
    msgs{end+1} = 'missing field sexclassifier_diagnostics_mean_nflies';
    success = false;
  else
    num_flies = round(metadata.sexclassifier_diagnostics_mean_nflies);
    if num_flies < check_params.min_num_flies,
      msgs{end+1} = sprintf('num_flies = round(%f) = %d < %d',...
        metadata.sexclassifier_diagnostics_mean_nflies,num_flies,check_params.min_num_flies);
      success = false;
    end
  end
    
  if ~isfield(metadata,'sexclassifier_diagnostics_mean_nfemales'),
    msgs{end+1} = 'missing field sexclassifier_diagnostics_mean_nfemales';
    success = false;
  else
    num_females = round(metadata.sexclassifier_diagnostics_mean_nfemales);
    if num_females < check_params.min_num_females,
      msgs{end+1} = sprintf('num_females = round(%f) = %d < %d',...
        metadata.sexclassifier_diagnostics_mean_nfemales,num_females,check_params.min_num_females);
      success = false;
    end
  end
    
  if ~isfield(metadata,'sexclassifier_diagnostics_mean_nmales'),
    msgs{end+1} = 'missing field sexclassifier_diagnostics_mean_nmales';
    success = false;
  else
    num_males = round(metadata.sexclassifier_diagnostics_mean_nmales);
    if num_males < check_params.min_num_males,
      msgs{end+1} = sprintf('num_males = round(%f) = %d < %d',...
        metadata.sexclassifier_diagnostics_mean_nmales,num_males,check_params.min_num_males);
      success = false;
    end
  end
  
end

%% check for Ctrax bugs: nan, inf in trajectories

ctraxfile = fullfile(expdir,dataloc_params.ctraxfilestr);
if ~exist(ctraxfile,'file'),
  msgs{end+1} = sprintf('Ctrax output mat file %s does not exist',ctraxfile);
  success = false;
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
      end
      
    end
  end
end

%% check for missing files


check_params.required_files = setdiff(check_params.required_files,...
  {'automatic_checks_incoming_results.txt'});
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
  end
end

%% output results to file, merging with automatedchecks incoming

if isempty(outfile_complete),
  fid = 1;
else
  fid = fopen(outfile_complete,'w');
end
if fid < 0,
  error('Could not open automatic checks results file %s for writing.',outfile);
end
if success && success_incoming,
  fprintf(fid,'automated_pf,P\n');
else
  fprintf(fid,'automated_pf,F\n');
  fprintf(fid,'notes_curation,');
  s = sprintf('%s\\n',msgs{:});
  s = s(1:end-2);
  if ~isempty(msgs_incoming),
    s_incoming = sprintf('%s\\n',msgs_incoming{:});
    s_incoming = s_incoming(1:end-2);
    if isempty(s),
      s = s_incoming;
    else
      s = [s_incoming,'\n',s];
    end
  end
  fprintf(fid,'%s\n',s);
end

if fid > 1,
  fclose(fid);
end

% for output

msgs = [msgs_incoming,msgs];
success = success && success_incoming;