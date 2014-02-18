function [isconsistent,filesmissing,msgs] = SAGEExpDirConsistencyCheck(expdir,varargin)

filesmissing = false;
isconsistent = true;
msgs = {};
missingfiles = {};

[dataset,analysis_protocol,settingsdir,datalocparamsfilestr,EPS,docheckhist,outfilename] = ...
  myparse(varargin,'dataset','data',...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'EPS',1e-3,...
  'docheckhist',false,...
  'outfilename','');

if expdir(end) == '/' || expdir(end) == '\',
  expdir = expdir(1:end-1);
end
[~,basename] = fileparts(expdir);
expname = ['FlyBowl_',basename];

%% get data from SAGE

data = SAGEGetBowlData('checkflags',false,'removemissingdata',false,...
  'dataset',dataset,'unflatten',false,...
  'experiment_name',expname,...
  'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol);

if docheckhist,
  data_hist = SAGEGetBowlData('checkflags',false,'removemissingdata',false,...
    'dataset','histogram','unflatten',false,...
    'experiment_name',expname,...
    'settingsdir',settingsdir,...
    'analysis_protocol',analysis_protocol);
end

if isempty(data),
  msgs{end+1} = sprintf('Experiment %s is not in SAGE',expname);  
  outfid = myfopen(outfilename,'a');
  fprintf(outfid,'\n%s\n',expdir);
  fprintf(outfid,'%s\n',msgs{:});
  myfclose(outfid);
  isconsistent = false;
  return;
end

%% get parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% check temperature stream

tempfile = fullfile(expdir,dataloc_params.temperaturefilestr);
if ~exist(tempfile,'file'),
  missingfiles{end+1} = sprintf('temperature stream file %s',dataloc_params.temperaturefilestr);
  filesmissing = true;
else
  tempdata = importdata(tempfile,',');
  if isempty(tempdata),
    msgs{end+1} = sprintf('Temperature stream data empty');
    filesmissing = true;
    missingfiles{end+1} = sprintf('temperature stream file %s',dataloc_params.temperaturefilestr);
  elseif size(tempdata,2) < 2,
    msgs{end+1} = sprintf('Temperature stream does not contain a second column');
    filesmissing = true;
    missingfiles{end+1} = sprintf('temperature stream file %s',dataloc_params.temperaturefilestr);
  else
    tempdata = tempdata(:,2);
    if ~isfield(data,'temperature_stream'),
      msgs{end+1} = sprintf('Temperature stream missing from SAGE');      
      isconsistent = false;
    elseif numel(tempdata) ~= numel(data.temperature_stream),
      msgs{end+1} = sprintf('Size mismatch for temperature stream');
      isconsistent = false;
    elseif max(abs(tempdata(:) - data.temperature_stream(:))) > EPS,
      msgs{end+1} = sprintf('Temperature stream data does not match');
      isconsistent = false;
    end
  end
end
      
%% check ufmf diagnostics

ufmffile = fullfile(expdir,dataloc_params.ufmfdiagnosticsfilestr);
if ~exist(ufmffile,'file'),
  %msgs{end+1} = sprintf('UFMF diagnostics file %s missing',ufmffile);
  missingfiles{end+1} = sprintf('UFMF diagnostics file %s',dataloc_params.ufmfdiagnosticsfilestr);
  filesmissing = true;
else
  [ufmfstats,success,errmsg] = readUFMFDiagnostics(ufmffile);
  if ~success,
    msgs{end+1} = sprintf('Error reading UFMF diagnostics from file: %s',errmsg);
    missingfiles{end+1} = sprintf('UFMF diagnostics file %s',dataloc_params.ufmfdiagnosticsfilestr);
    filesmissing = true;
  else
    % check summary
    
    % some of the field names are slightly different
    dict = {'fracFramesWithFracFgPx_gt_0_050000','fracFramesWithFracFgPx0050000'
      'fracFramesWithFracFgPx_gt_0_100000','fracFramesWithFracFgPx0100000'
      'fracFramesWithFracFgPx_gt_0_250000','fracFramesWithFracFgPx0250000'};
    [isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(ufmfstats.summary,'ufmf_diagnostics_summary',data,dict);
    isconsistent = isconsistent && isconsistent1;
    filesmissing = filesmissing || filesmissing1;
    msgs = [msgs,msgs1];
    
    % check stream
    [isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(ufmfstats.stream,'ufmf_diagnostics_stream',data,dict);
    isconsistent = isconsistent && isconsistent1;
    filesmissing = filesmissing || filesmissing1;
    msgs = [msgs,msgs1];

  end
    
end

%% check quickstats

quickfile = fullfile(expdir,dataloc_params.quickstatsfilestr);
[isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(quickfile,'QuickStats',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;
msgs = [msgs,msgs1];
if filesmissing1,
  missingfiles{end+1} = sprintf('QuickStats file %s',dataloc_params.quickstatsfilestr);
end

%% check automated_checks_incoming_results

incomingfile = fullfile(expdir,dataloc_params.automaticchecksincomingresultsfilestr);
success = false;
if ~exist(incomingfile,'file'),
  %msgs{end+1} = sprintf('Automatic checks incoming results file %s does not exist',incomingfile);
  filesmissing = true;
  missingfiles{end+1} = sprintf('Automatic checks incoming results file %s',dataloc_params.automaticchecksincomingresultsfilestr);

else
  try
    incoming_results = ReadParams(incomingfile);
    success = true;
  catch ME,
    msgs{end+1} = sprintf('Could not read Ctrax diagnostics file %s',incomingfile);
    getReport(ME)
    filesmissing = true;
    missingfiles{end+1} = sprintf('Automatic checks incoming results file %s',dataloc_params.automaticchecksincomingresultsfilestr);
  end
end

if success,
  if ~isfield(incoming_results,'automated_pf'),
    msgs{end+1} = sprintf('automated_pf not set in incoming checks results file %s',incomingfile);
  elseif ~isfield(data,'automated_pf'),
    msgs{end+1} = sprintf('automated_pf missing from SAGE');
    isconsistent = false;
  elseif (incoming_results.automated_pf == 'F') && (data.automated_pf ~= 'F'),
    msgs{end+1} = sprintf('automated_pf = %s in SAGE mismatches %s in %s',data.automated_pf,incoming_results.automated_pf,incomingfile);
    isconsistent = false;
  end
  
  if isfield(incoming_results,'notes_curation'),
    if ~isfield(data,'notes_curation'),
      msgs{end+1} = sprintf('notes_curation missing from SAGE');
      isconsistent = false;
    else
      expr = ['.*',strrep(incoming_results.notes_curation,'\n','(\n)*\s*'),'.*'];
      if isempty(regexp(data.notes_curation,expr,'once')),
        msgs{end+1} = sprintf('notes_curation string %s in %s, but missing from SAGE''s notes_curation = %s',incoming_results.notes_curation,incomingfile,data.notes_curation);
        isconsistent = false;
      end
    end    
  end
end

%% check ctrax_diagnostics

ctraxfile = fullfile(expdir,dataloc_params.ctraxdiagnosticsfilestr);
[isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(ctraxfile,'ctrax_diagnostics',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;
msgs = [msgs,msgs1];
if filesmissing1,
  missingfiles{end+1} = sprintf('ctrax diagnostics file %s',dataloc_params.ctraxdiagnosticsfilestr);
end

%% check registration data

regfile = fullfile(expdir,dataloc_params.registrationtxtfilestr);
[isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(regfile,'registrationdata',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;
msgs = [msgs,msgs1];
if filesmissing1,
  missingfiles{end+1} = sprintf('registration file %s',dataloc_params.registrationtxtfilestr);
end

%% check sexclassifier data

sexfile = fullfile(expdir,dataloc_params.sexclassifierdiagnosticsfilestr);
[isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(sexfile,'sexclassifier_diagnostics',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;
msgs = [msgs,msgs1];
if filesmissing1,
  missingfiles{end+1} = sprintf('sex classifier diagnostics file %s',dataloc_params.sexclassifierdiagnosticsfilestr);
end

%% check perframestats

statsfile = fullfile(expdir,dataloc_params.statsperframetxtfilestr);
[isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(statsfile,'stats',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;
msgs = [msgs,msgs1];
if filesmissing1,
  missingfiles{end+1} = sprintf('stats perframe text file %s',dataloc_params.statsperframetxtfilestr);
end

%% check histogram data

if docheckhist,
  histfile = fullfile(expdir,dataloc_params.histperframetxtfilestr);
  [isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(histfile,'hist',data_hist);
  isconsistent = isconsistent && isconsistent1;
  filesmissing = filesmissing || filesmissing1;
  msgs = [msgs,msgs1];
  if filesmissing1,
    missingfiles{end+1} = sprintf('histogram per frame file %s',dataloc_params.histperframetxtfilestr);
  end
end
  
%% check bias diagnostics

biasfile = fullfile(expdir,dataloc_params.biasdiagnosticsfilestr);
[isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(biasfile,'bias_diagnostics',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;
msgs = [msgs,msgs1];
if filesmissing1,
  missingfiles{end+1} = sprintf('bias diagnostics file %s',dataloc_params.biasdiagnosticsfilestr);
end

%% check bkgd diagnostics

bkgdfile = fullfile(expdir,dataloc_params.bkgddiagnosticsfilestr);
[isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(bkgdfile,'bkgd_diagnostics',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;
msgs = [msgs,msgs1];
if filesmissing1,
  missingfiles{end+1} = sprintf('bkgd diagnostics file %s',dataloc_params.bkgddiagnosticsfilestr);
end

%% check temperature diagnostics

temperaturefile = fullfile(expdir,dataloc_params.temperaturediagnosticsfilestr);
[isconsistent1,filesmissing1,msgs1] = CheckDiagnosticsFileSAGEConsistency(temperaturefile,'temperature_diagnostics',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;
msgs = [msgs,msgs1];
if filesmissing1,
  missingfiles{end+1} = sprintf('temperature diagnostics file %s',dataloc_params.temperaturediagnosticsfilestr);
end

%% output to file

if ~isconsistent,
  outfid = myfopen(outfilename,'a');
  fprintf(outfid,'%s\t',expdir);
  fprintf(outfid,'%s, ',msgs{:});
  if filesmissing,
    fprintf(outfid,'\tMissing files: ');
    fprintf(outfid,'%s, ',missingfiles{:});
  else
    fprintf(outfid,'\t');
  end
  fnsprint = {'automated_pf','automated_pf_category','manual_pf','screen_type'};
  for i = 1:numel(fnsprint),
    fn = fnsprint{i};
    if isfield(data,fn),
      fprintf(outfid,'%s = %s\t',fn,data.(fn));
    else
      fprintf(outfid,'missing %s\t',fn);
    end
  end  
  fprintf(outfid,'\n');
  myfclose(outfid);
elseif filesmissing && incoming_results.automated_pf ~= 'F',
  fprintf('\n%s consistent but missing files:\n',expdir);
  fprintf('%s\n',missingfiles{:});
end