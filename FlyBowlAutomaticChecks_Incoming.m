function [success,msgs,iserror] = FlyBowlAutomaticChecks_Incoming(expdir,varargin)

success = true;
msgs = {};

datetime_format = 'yyyymmddTHHMMSS';

[analysis_protocol,settingsdir,datalocparamsfilestr,DEBUG,min_barcode_expdatestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt','debug',false,...
  'min_barcode_expdatestr','20110301T000000');
min_barcode_expdatenum = datenum(min_barcode_expdatestr,datetime_format);

%% parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
paramsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.automaticchecksincomingparamsfilestr);
check_params = ReadParams(paramsfile);
if ~iscell(check_params.control_line_names),
  check_params.control_line_names = {check_params.control_line_names};
end
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
ufmfdiagnosticsfile = fullfile(expdir,dataloc_params.ufmfdiagnosticsfilestr);
temperaturefile = fullfile(expdir,dataloc_params.temperaturefilestr);
outfile = fullfile(expdir,dataloc_params.automaticchecksincomingresultsfilestr);

% order matters here: higher up categories have higher priority
categories = {'flag_aborted_set_to_1',...
  'missing_video',...
  'missing_metadata_fields',...
  'short_video',...
  'missing_capture_files',...
  'flag_redo_set_to_1',...
  'fliesloaded_time_too_short',...
  'fliesloaded_time_too_long',...
  'shiftflytemp_time_too_long',...
  'no_barcode',...
  'incoming_checks_other'};
category2idx = struct;
for i = 1:numel(categories),
  category2idx.(categories{i}) = i;
end
iserror = false(1,numel(categories));

%% read metadata

metadata = ReadMetadataFile(metadatafile);

%% check for metadata fields

required_fns = {'flag_aborted','flag_redo','seconds_fliesloaded','seconds_shiftflytemp',...
  'screen_type','line','cross_barcode'};
ismissingfn = ~ismember(required_fns,fieldnames(metadata));
if any(ismissingfn),
  success = false;
  msgs{end+1} = ['Missing metadata fields:',sprintf(' %s',required_fns{ismissingfn})];
  iserror(category2idx.missing_metadata_fields) = true;
end

%% check for flags

if metadata.flag_aborted ~= 0,
  success = false;
  msgs{end+1} = 'Experiment aborted.';
  iserror(category2idx.flag_aborted_set_to_1) = true;
end

if metadata.flag_redo ~= 0,
  success = false;
  msgs{end+1} = 'Redo flag set to 1.';
  iserror(category2idx.flag_redo_set_to_1) = true;
end

%% check loading time

if metadata.seconds_fliesloaded < check_params.min_seconds_fliesloaded,
  success = false;
  msgs{end+1} = sprintf('Load time = %f < %f seconds.',metadata.seconds_fliesloaded,check_params.min_seconds_fliesloaded);
  iserror(category2idx.fliesloaded_time_too_short) = true;
end
if metadata.seconds_fliesloaded > check_params.max_seconds_fliesloaded,
  success = false;
  msgs{end+1} = sprintf('Load time = %f > %f seconds.',metadata.seconds_fliesloaded,check_params.max_seconds_fliesloaded);
  iserror(category2idx.fliesloaded_time_too_long) = true;
end

%% check shiftflytemp time
if metadata.seconds_shiftflytemp > check_params.max_seconds_shiftflytemp,
  success = false;
  msgs{end+1} = sprintf('Shift fly temp time = %f > %f seconds.',metadata.seconds_shiftflytemp,check_params.max_seconds_shiftflytemp);
  iserror(category2idx.shiftflytemp_time_too_long) = true;
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

exp_datenum = datenum(metadata.exp_datetime,datetime_format);
  
if (~ismember(metadata.line,check_params.control_line_names) || ...
    (exp_datenum >= min_barcode_expdatenum)) && ...
    metadata.cross_barcode < 0,
  success = false;
  msgs{end+1} = 'Barcode = -1 and line_name indicates not control';
  iserror(category2idx.no_barcode) = true;
end

%% check video length

if ~exist(ufmfdiagnosticsfile,'file'),
  success = false;
  msgs{end+1} = 'Missing UFMF diagnostics file.';
  iserror(category2idx.missing_capture_files) = true;
else
  [ufmf_diagnostics,success1,errmsg1] = readUFMFDiagnostics(ufmfdiagnosticsfile);
  if ~success1,
    success = false;
    msgs{end+1} = sprintf('Error reading UFMF diagnostics file: %s',errmsg1);
    iserror(category2idx.incoming_checks_other) = true;
  else
    if ~isfield(ufmf_diagnostics,'summary') || ~isfield(ufmf_diagnostics.summary,'nFrames'),
      success = false;
      msgs{end+1} = 'ufmf_diagnostics_summary_nFrames missing.';
      iserror(category2idx.incoming_checks_other) = true;
    else
      if ufmf_diagnostics.summary.nFrames < check_params.min_ufmf_diagnostics_summary_nframes,
        success = false;
        msgs{end+1} = sprintf('Video contains %d < %d frames.',ufmf_diagnostics.summary.nFrames,check_params.min_ufmf_diagnostics_summary_nframes);
        iserror(category2idx.short_video) = true;
      end
    end
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

fn = dataloc_params.moviefilestr;
isfile = exist(fullfile(expdir,fn),'file');
if ~isfile,
  iserror(category2idx.missing_video) = true;
end

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
    iserror(category2idx.missing_capture_files) = true;
  end
end

%% output results to file

if DEBUG,
  fid = 1;
else
  fid = fopen(outfile,'w');
end
if fid < 0,
  error('Could not open automatic checks results file %s for writing.',outfile);
end
if success,
  fprintf(fid,'automated_pf,U\n');
else
  fprintf(fid,'automated_pf,F\n');
  fprintf(fid,'notes_curation,');
  s = sprintf('%s\\n',msgs{:});
  s = s(1:end-2);
  fprintf(fid,'%s\n',s);
  i = find(iserror,1);
  if isempty(i),
    s = 'incoming_checks_other';
  else
    s = categories{i};
  end
  fprintf(fid,'automated_pf_category,%s\n',s);
end

if ~DEBUG,
  fclose(fid);
end