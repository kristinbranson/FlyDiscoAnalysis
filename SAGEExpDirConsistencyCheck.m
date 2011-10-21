function [isconsistent,filesmissing] = SAGEExpDirConsistencyCheck(expdir,varargin)

filesmissing = false;
isconsistent = true;

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

if isempty(outfilename),
  outfid = 1;
else
  outfid = fopen(outfilename,'a');
  if outfid < 0,
    error('Could not open file %s for appending\n',outfilename);
  end
end

fprintf(outfid,'\n%s\n',expdir);

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
  fprintf(outfid,'Experiment %s is not in SAGE\n',expname);
  isconsistent = false;
  return;
end

%% get parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% check temperature stream

tempfile = fullfile(expdir,dataloc_params.temperaturefilestr);
if ~exist(tempfile,'file'),
  fprintf(outfid,'Temperature stream file %s missing\n',tempfile);
  filesmissing = true;
else
  tempdata = importdata(tempfile,',');
  if isempty(tempdata),
    fprintf(outfid,'Temperature stream data empty\n');
  elseif size(tempdata,2) < 2,
    fprintf(outfid,'Temperature stream does not contain a second column\n');
  else
    tempdata = tempdata(:,2);
    if ~isfield(data,'temperature_stream'),
      fprintf(outfid,'Temperature stream missing from SAGE\n');      
      isconsistent = false;
    elseif numel(tempdata) ~= numel(data.temperature_stream),
      fprintf(outfid,'Size mismatch for temperature stream\n');
      isconsistent = false;
    elseif max(abs(tempdata(:) - data.temperature_stream(:))) > EPS,
      fprintf(outfid,'Temperature stream data does not match\n');
      isconsistent = false;
    end
  end
end
      
%% check ufmf diagnostics

ufmffile = fullfile(expdir,dataloc_params.ufmfdiagnosticsfilestr);
if ~exist(ufmffile,'file'),
  fprintf(outfid,'UFMF diagnostics file %s missing\n',ufmffile);
  filesmissing = true;
else
  [ufmfstats,success,errmsg] = readUFMFDiagnostics(ufmffile);
  if ~success,
    fprintf(outfid,'Error reading UFMF diagnostics from file: %s\n',errmsg);
  else
    % check summary
    
    % some of the field names are slightly different
    dict = {'fracFramesWithFracFgPx_gt_0_050000','fracFramesWithFracFgPx0050000'
      'fracFramesWithFracFgPx_gt_0_100000','fracFramesWithFracFgPx0100000'
      'fracFramesWithFracFgPx_gt_0_250000','fracFramesWithFracFgPx0250000'};
    [isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,ufmfstats.summary,'ufmf_diagnostics_summary',data,dict);
    isconsistent = isconsistent && isconsistent1;
    filesmissing = filesmissing || filesmissing1;
    
    % check stream
    [isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,ufmfstats.stream,'ufmf_diagnostics_stream',data,dict);
    isconsistent = isconsistent && isconsistent1;
    filesmissing = filesmissing || filesmissing1;

  end
    
end

%% check quickstats

quickfile = fullfile(expdir,dataloc_params.quickstatsfilestr);
[isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,quickfile,'QuickStats',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;


%% check automated_checks_incoming_results

incomingfile = fullfile(expdir,dataloc_params.automaticchecksincomingresultsfilestr);
success = false;
if ~exist(incomingfile,'file'),
  fprintf(outfid,'Automatic checks incoming results file %s does not exist\n',incomingfile);
  filesmissing = true;
else
  try
    incoming_results = ReadParams(incomingfile);
    success = true;
  catch ME,
    fprintf(outfid,'Could not read Ctrax diagnostics file %s\n',incomingfile);
    getReport(ME)
    filesmissing = true;
  end
end

if success,
  if ~isfield(incoming_results,'automated_pf'),
    fprintf(outfid,'automated_pf not set in incoming checks results file %s\n',incomingfile);
  elseif ~isfield(data,'automated_pf'),
    fprintf(outfid,'automated_pf missing from SAGE\n');
    isconsistent = false;
  elseif (incoming_results.automated_pf == 'F') && (data.automated_pf ~= 'F'),
    fprintf(outfid,'automated_pf = %s in SAGE mismatches %s in %s\n',data.automated_pf,incoming_results.automated_pf,incomingfile);
    isconsistent = false;
  end
  
  if isfield(incoming_results,'notes_curation'),
    if ~isfield(data,'notes_curation'),
      fprintf(outfid,'notes_curation missing from SAGE\n');
      isconsistent = false;
    else
      expr = ['.*',strrep(incoming_results.notes_curation,'\n','(\n)*\s*'),'.*'];
      if isempty(regexp(data.notes_curation,expr,'once')),
        fprintf(outfid,'notes_curation string %s in %s, but missing from SAGE''s notes_curation = %s\n',incoming_results.notes_curation,incomingfile,data.notes_curation);
        isconsistent = false;
      end
    end    
  end
end

%% check ctrax_diagnostics

ctraxfile = fullfile(expdir,dataloc_params.ctraxdiagnosticsfilestr);
[isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,ctraxfile,'ctrax_diagnostics',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;

%% check registration data

regfile = fullfile(expdir,dataloc_params.registrationtxtfilestr);
[isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,regfile,'registrationdata',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;

%% check sexclassifier data

sexfile = fullfile(expdir,dataloc_params.sexclassifierdiagnosticsfilestr);
[isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,sexfile,'sexclassifier_diagnostics',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;

%% check perframestats

statsfile = fullfile(expdir,dataloc_params.statsperframetxtfilestr);
[isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,statsfile,'stats',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;

%% check histogram data

if docheckhist,
  histfile = fullfile(expdir,dataloc_params.histperframetxtfilestr);
  [isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,histfile,'hist',data_hist);
  isconsistent = isconsistent && isconsistent1;
  filesmissing = filesmissing || filesmissing1;
end
  
%% check bias diagnostics

biasfile = fullfile(expdir,dataloc_params.biasdiagnosticsfilestr);
[isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,biasfile,'bias_diagnostics',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;

%% check bkgd diagnostics

bkgdfile = fullfile(expdir,dataloc_params.bkgddiagnosticsfilestr);
[isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,bkgdfile,'bkgd_diagnostics',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;

%% check temperature diagnostics

temperaturefile = fullfile(expdir,dataloc_params.temperaturediagnosticsfilestr);
[isconsistent1,filesmissing1] = CheckDiagnosticsFileSAGEConsistency(outfid,temperaturefile,'temperature_diagnostics',data);
isconsistent = isconsistent && isconsistent1;
filesmissing = filesmissing || filesmissing1;

%% close output file

if ~isempty(outfilename),
  fclose(outfid);
end