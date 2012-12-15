function [success,msgs,stage] = FlyBowlAnalysisPipeline_heberlein(expdir,varargin)

version = '0.1';
timestamp = datestr(now,'yyyymmddTHHMMSS');

% expdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_housing_pBDPGAL4U_20111216/results/pBDPGAL4U_UAS_IVS_myr_3_0009_Rig1Plate15BowlA_20120113T142907';
% expdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_mating_galit_CS_20120211/results/EXT_CantonS_1220002_None_Rig1Plate15BowlA_20120211T123840';
success = false;
msgs = {};
stage = 'start';

[analysis_protocol,settingsdir,datalocparamsfilestr,...
  automaticchecksincoming_params,registration_params,sexclassification_params,...
  computeperframefeatures_params,computeperframestats_params,...
  plotperframestats_params,makectraxresultsmovie_params,...
  extradiagnostics_params,automaticcheckscomplete_params,...
  forcecompute,doautomaticchecksincoming,doregistration,dosexclassification,...
  dotrackwings,docomputeperframefeatures,dojaabadetect,docomputeperframestats,doplotperframestats,...
  domakectraxresultsmovie,doextradiagnostics,doanalysisprotocol,...
  doautomaticcheckscomplete,...
  requiredfiles_start,requiredfiles_automaticchecksincoming,...
  requiredfiles_registration,...
  requiredfiles_sexclassification,...
  requiredfiles_wingtracking,...
  requiredfiles_computeperframefeatures,...
  requiredfiles_computeperframestats,...
  requiredfiles_plotperframestats,...
  requiredfiles_makectraxresultsmovie,...
  requiredfiles_extradiagnostics,...
  requiredfiles_analysisprotocol,...
  requiredfiles_automaticcheckscomplete] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'automaticchecksincoming_params',{},...
  'registration_params',{},...
  'sexclassification_params',{},...
  'computeperframefeatures_params',{},...
  'computeperframestats_params',{},...
  'plotperframestats_params',{},...
  'makectraxresultsmovie_params',{},...
  'extradiagnostics_params',{},...
  'automaticcheckscomplete_params',{},...
  'forcecompute',false,...
  'doautomaticchecksincoming',true,...
  'doregistration',true,...
  'dosexclassification',true,...
  'dotrackwings',false,...
  'docomputeperframefeatures',true,...
  'dojaabadetect',false,...
  'docomputeperframestats',false,...
  'doplotperframestats',false,...
  'domakectraxresultsmovie',true,...
  'doextradiagnostics',true,...
  'doanalysisprotocol',isunix,...
  'doautomaticcheckscomplete',true,...
  'requiredfiles_start',{'annfilestr','ctraxfilestr'},...
  'requiredfiles_automaticchecksincoming',{'automaticchecksincomingresultsfilestr'},...
  'requiredfiles_registration',{'trxfilestr','registrationmatfilestr','registrationtxtfilestr','registrationimagefilestr'},...
  'requiredfiles_sexclassification',{'sexclassifierdiagnosticsfilestr'},...
  'requiredfiles_wingtracking',{'sexclassifierdiagnosticsfilestr'},...
  'requiredfiles_computeperframefeatures',{'perframedir','PERFRAMEMATFILES'},...
  'requiredfiles_computeperframestats',{'statsperframetxtfilestr','statsperframematfilestr','histperframetxtfilestr','histperframematfilestr'},...
  'requiredfiles_plotperframestats',{'figdir','ANALYSISPLOTS'},...
  'requiredfiles_makectraxresultsmovie',{'CTRAXRESULTSMOVIE'},...
  'requiredfiles_extradiagnostics',{'biasdiagnosticsimagefilestr','biasdiagnosticsfilestr','biasdiagnosticsmatfilestr',...
  'temperaturediagnosticsfilestr','bkgddiagnosticsimagefilestr','bkgddiagnosticsfilestr','bkgddiagnosticsmatfilestr',...
  'videodiagnosticsfilestr','videodiagnosticsmatfilestr','videodiagnosticsimagefilestr'},...
  'requiredfiles_analysisprotocol',{'analysis_protocol.txt'},...
  'requiredfiles_automaticcheckscomplete',{'automaticcheckscompleteresultsfilestr'});

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

% open log file
if isfield(dataloc_params,'pipeline_logfilestr'),
  logfile = fullfile(expdir,dataloc_params.pipeline_logfilestr);
  logfid = fopen(logfile,'a');
  if logfid < 1,
    warning('Could not open log file %s\n',logfile);
    logfid = 1;
  end
else
  logfid = 1;
end

fprintf(logfid,'\n\n***\nRunning FlyBowlAnalysisPipeline_heberlein version %s analysis_protocol %s at %s\n',version,analysis_protocol,timestamp);

if ischar(forcecompute),
  forcecompute = str2double(forcecompute) ~= 0;
end
fprintf(logfid,'forcecompute: %d\n',forcecompute);

if ischar(doautomaticchecksincoming),
  doautomaticchecksincoming = str2double(doautomaticchecksincoming) ~= 0;
end
fprintf(logfid,'doautomaticchecksincoming: %d\n',doautomaticchecksincoming);

if ischar(doregistration),
  doregistration = str2double(doregistration) ~= 0;
end
fprintf(logfid,'doregistration: %d\n',doregistration);

if ischar(dosexclassification),
  dosexclassification = str2double(dosexclassification) ~= 0;
end
fprintf(logfid,'dosexclassification: %d\n',dosexclassification);

if ischar(docomputeperframefeatures),
  docomputeperframefeatures = str2double(docomputeperframefeatures) ~= 0;
end
fprintf(logfid,'docomputeperframefeatures: %d\n',docomputeperframefeatures);

if ischar(dojaabadetect),
  dojaabadetect = str2double(dojaabadetect) ~= 0;
end
fprintf(logfid,'dojaabadetect: %d\n',dojaabadetect);

if ischar(docomputeperframestats),
  docomputeperframestats = str2double(docomputeperframestats) ~= 0;
end
fprintf(logfid,'docomputeperframestats: %d\n',docomputeperframestats);

if ischar(doplotperframestats),
  doplotperframestats = str2double(doplotperframestats) ~= 0;
end
fprintf(logfid,'doplotperframestats: %d\n',doplotperframestats);

if ischar(domakectraxresultsmovie),
  domakectraxresultsmovie = str2double(domakectraxresultsmovie) ~= 0;
end
fprintf(logfid,'domakectraxresultsmovie: %d\n',domakectraxresultsmovie);

if ischar(doextradiagnostics),
  doextradiagnostics = str2double(doextradiagnostics) ~= 0;
end
fprintf(logfid,'doextradiagnostics: %d\n',doextradiagnostics);

if ischar(doautomaticcheckscomplete),
  doautomaticcheckscomplete = str2double(doautomaticcheckscomplete) ~= 0;
end
fprintf(logfid,'doautomaticcheckscomplete: %d\n',doautomaticcheckscomplete);

%% check that experiment exists

if ~exist(expdir,'dir'),
  msgs = {sprintf('Experiment directory %s does not exist',expdir)};
  fprintf(logfid,'%s\n',msgs{:});
  if logfid > 1, fclose(logfid); end
  return;
end

%% check for files required at start

[ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_start);
if ismissingfile,
  msgs = cellfun(@(x) sprintf('Missing start file %s',x),missingfiles,'UniformOutput',false);
  fprintf(logfid,'%s\n',msgs{:});
  if logfid > 1, fclose(logfid); end
  return;
end

%% incoming checks

stage = 'automaticchecks_incoming';

if doautomaticchecksincoming,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticchecksincoming);
  if forcecompute || todo,
    try
      fprintf(logfid,'AutomaticChecks_Incoming...\n');

      [success1,msgs] = FlyBowlAutomaticChecks_Incoming(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        automaticchecksincoming_params{:});
      
      if ~success1,
        fprintf(logfid,'AutomaticChecks_Incoming failed:\n');
        fprintf(logfid,'%s\n',msgs{:});
        if logfid > 1, fclose(logfid); end
        return;
      end
    catch ME,
      msgs = {sprintf('Error running AutomaticChecks_Incoming:\n%s',getReport(ME))};
      fprintf(logfid,'AutomaticChecks_Incoming failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end
      return;
    end
  end
  
  % make sure automatic checks files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticchecksincoming);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing automaticchecks_incoming file %s',x),missingfiles,'UniformOutput',false);
    fprintf(logfid,'AutomaticChecks_Incoming failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    if logfid > 1, fclose(logfid); end; return;
  end
  
end

%% registration

stage = 'registration';

if doregistration,
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_registration);
  if forcecompute || todo,
    
    try
      fprintf(logfid,'RegisterTrx...\n');

      FlyBowlRegisterTrx(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        registration_params{:});
      
    catch ME,
      msgs = {sprintf('Error running RegisterTrx:\n%s',getReport(ME))};
      fprintf(logfid,'RegisterTrx failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end; return;
    end
  end
  
  % make sure registration files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_registration);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing registration file %s',x),missingfiles,'UniformOutput',false);
    fprintf(logfid,'RegisterTrx failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    if logfid > 1, fclose(logfid); end; return;
  end
  
end

%% sex classification


stage = 'sexclassification';

if dosexclassification,
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_sexclassification);
  if forcecompute || todo,
    
    try
      fprintf(logfid,'ClassifySex...\n');
      FlyBowlClassifySex2(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        sexclassification_params{:});
      
    catch ME,
      msgs = {sprintf('Error running SexClassification:\n%s',getReport(ME))};
      fprintf(logfid,'SexClassification failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end; return;
    end
  end
  
  % make sure sexclassification files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_sexclassification);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing sexclassification file %s',x),missingfiles,'UniformOutput',false);
    fprintf(logfid,'SexClassification failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    if logfid > 1, fclose(logfid); end; return;
  end
  
end

%% wing tracking

stage = 'trackwings';
if dotrackwings,
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_wingtracking);
  if forcecompute || todo,
    
    try
      fprintf(logfid,'Wing tracking...\n');
      TrackWings(expdir,...
        'moviefilestr',dataloc_params.moviefilestr,...
        'annfilestr',dataloc_params.annfilestr,...
        'trxfilestr',dataloc_params.trxfilestr,...
        'outtrxfilestr',dataloc_params.wingtrxfilestr,...
        'perframedir',dataloc_params.perframedir,...
        'paramsfile',fullfile(settingsdir,analysis_protocol,dataloc_params.wingtrackingparamsfile));
      
    catch ME,
      msgs = {sprintf('Error running TrackWings:\n%s',getReport(ME))};
      fprintf(logfid,'TrackWings failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end; return;
    end
  end
  
end

%% compute per-frame features

stage = 'computeperframefeatures';

if docomputeperframefeatures,
    
  if ismember('PERFRAMEMATFILES',requiredfiles_computeperframefeatures),
    % read in per-frame fns to compute
    i = find(strcmpi('perframefns',computeperframefeatures_params),1);
    if ~isempty(i),
      perframefns = computeperframefeatures_params{i+1};
    else
      perframefnsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.perframefnsfilestr);
      perframefns = importdata(perframefnsfile);      
    end
    perframematfiles = cellfun(@(x) fullfile(dataloc_params.perframedir,[x,'.mat']),perframefns,'UniformOutput',false);
    requiredfiles_computeperframefeatures = [requiredfiles_computeperframefeatures(:)',perframematfiles(:)'];
    requiredfiles_computeperframefeatures = setdiff(requiredfiles_computeperframefeatures,{'PERFRAMEMATFILES'});
  end
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframefeatures);
  if forcecompute || todo,
    
    try
      fprintf(logfid,'ComputePerFrameFeatures...\n');
      FlyBowlComputePerFrameFeatures(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        'forcecompute',forcecompute,...
        computeperframefeatures_params{:});
      
    catch ME,
      msgs = {sprintf('Error running ComputePerFrameFeatures:\n%s',getReport(ME))};
      fprintf(logfid,'ComputePerFrameFeatures failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end; return;
    end
  end
  
  % make sure computeperframefeatures files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframefeatures);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing computeperframefeatures file %s',x),missingfiles,'UniformOutput',false);
    fprintf(logfid,'ComputePerFrameFeatures failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    if logfid > 1, fclose(logfid); end; return;
  end
  
end

%% behavior detection

stage = 'jaabadetect';
if dojaabadetect,
  
  try

  pwdprev = pwd;  

    
  classifierparamsfiles = dataloc_params.classifierparamsfiles;
  classifierparams = cell(1,numel(classifierparamsfiles));
  jaabadir = fileparts(which('JAABADetect'));
  cd(jaabadir);
  for i = 1:numel(classifierparamsfiles),
    classifierparamsfile = fullfile(settingsdir,analysis_protocol,classifierparamsfiles{i});
    classifierparams{i} = ReadClassifierParamsFile(classifierparamsfile);
    % remove behaviors that are already done
    if ~forcecompute,
      issuccess = CheckScores({expdir},classifierparamsfiles,classifierparams{i});
      classifierparams{i}(issuccess) = [];
    end
  end
   
  for i = 1:numel(classifierparams),
    if isempty(classifierparams{i}),
      continue;
    end
    classifierfiles = {classifierparams{i}.classifierfile};
    configparams = rmfield(classifierparams{i},{'classifierfile','configfile'});
    JAABADetect(expdir,'classifierfiles',classifierfiles,'configparams',configparams);
  end
  
  cd(pwdprev);
  
  catch ME,
    msgs = {sprintf('Error running JAABADetect:\n%s',getReport(ME))};
    fprintf(logfid,'JAABADetect failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    cd(pwdprev);
    if logfid > 1, fclose(logfid); end; return;
  end
  
end

%% compute per-frame statistics

stage = 'computeperframestats';

if docomputeperframestats,
    
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframestats);
  if forcecompute || todo,
    
    try
      fprintf(logfid,'ComputePerFrameStats2...\n');
      FlyBowlComputePerFrameStats2(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        computeperframestats_params{:});
      
    catch ME,
      msgs = {sprintf('Error running ComputePerFrameStats:\n%s',getReport(ME))};
      fprintf(logfid,'ComputePerFrameStats failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end; return;
    end
  end
  
  % make sure computeperframestats files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframestats);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing computeperframestats file %s',x),missingfiles,'UniformOutput',false);
    fprintf(logfid,'ComputePerFrameStats failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    if logfid > 1, fclose(logfid); end; return;
  end
  
end

%% plot per-frame stats

stage = 'plotperframestats';

if doplotperframestats,
    
  if ismember('ANALYSISPLOTS',requiredfiles_plotperframestats),
    histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
    hist_perframefeatures = ReadHistPerFrameFeatures2(histperframefeaturesfile);
    files = cell(1,numel(hist_perframefeatures));
    for i = 1:numel(hist_perframefeatures),
      if strcmp(hist_perframefeatures(i).framecondition,'any'),
        files{i} = fullfile(dataloc_params.figdir,sprintf('hist_%s.png',hist_perframefeatures(i).field));
      else
        files{i} = fullfile(dataloc_params.figdir,sprintf('hist_%s_%s.png',hist_perframefeatures(i).field,hist_perframefeatures(i).framecondition));
      end
    end
    requiredfiles_plotperframestats = unique(files);
%     
%     hist_fields = unique({hist_perframefeatures.field});
%     for i = 1:numel(hist_fields),
%       savename = sprintf('hist_%s.png',hist_fields{i});
%       savename = fullfile(dataloc_params.figdir,savename);
%       requiredfiles_plotperframestats{end+1} = savename; %#ok<AGROW>
%     end
    requiredfiles_plotperframestats = setdiff(requiredfiles_plotperframestats,{'ANALYSISPLOTS'});
    requiredfiles_plotperframestats{end+1} = fullfile(dataloc_params.figdir,'stats.png');
  end
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_plotperframestats);
  if forcecompute || todo,
    
    try
      fprintf(logfid,'PlotPerFrameStats...\n');
      FlyBowlPlotPerFrameStats2(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        plotperframestats_params{:});
    catch ME,
      msgs = {sprintf('Error running PlotPerFrameStats:\n%s',getReport(ME))};
      fprintf(logfid,'PlotPerFrameStats failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end; return;
    end
  end
  
  % make sure plotperframestats files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_plotperframestats);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing plotperframestats file %s',x),missingfiles,'UniformOutput',false);
    fprintf(logfid,'PlotPerFrameStats failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    if logfid > 1, fclose(logfid); end; return;
  end
  
end

%% make results movie

stage = 'ctraxresultsmovie';

if domakectraxresultsmovie,
    
  i = find(strcmp('CTRAXRESULTSMOVIE',requiredfiles_makectraxresultsmovie),1);
  if ~isempty(i),
    [~,basename] = fileparts(expdir);
    avifilestr = sprintf('%s_%s',dataloc_params.ctraxresultsavifilestr,basename);
    xvidfile = [avifilestr,'.avi'];
    requiredfiles_makectraxresultsmovie{i} = xvidfile;
  end
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makectraxresultsmovie);
  if forcecompute || todo,
    
    try
      fprintf(logfid,'MakeCtraxResultsMovie...\n');
      FlyBowlMakeCtraxResultsMovie(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        makectraxresultsmovie_params{:});
    catch ME,
      msgs = {sprintf('Error running MakeCtraxResultsMovie:\n%s',getReport(ME))};
      fprintf(logfid,'MakeCtraxResultsMovie failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end; return;
    end
  end
  
  % make sure makectraxresultsmovie files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makectraxresultsmovie);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing makectraxresultsmovie file %s',x),missingfiles,'UniformOutput',false);
    fprintf(logfid,'MakeCtraxResultsMovie failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    if logfid > 1, fclose(logfid); end; return;
  end
  
end

%% extra diagnostics

stage = 'extradiagnostics';

if doextradiagnostics,
    
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_extradiagnostics);
  if forcecompute || todo,
    
    try
      fprintf(logfid,'ExtraDiagnostics...\n');
      FlyBowlExtraDiagnostics(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        extradiagnostics_params{:});
    catch ME,
      msgs = {sprintf('Error running ExtraDiagnostics:\n%s',getReport(ME))};
      fprintf(logfid,'ExtraDiagnostics failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end; return;
    end
  end
  
  % make sure extradiagnostics files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_extradiagnostics);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing extradiagnostics file %s',x),missingfiles,'UniformOutput',false);
    fprintf(logfid,'ExtraDiagnostics failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    if logfid > 1, fclose(logfid); end; return;
  end
  
end

%% generate analysis protocol file

stage = 'analysisprotocol';

if doanalysisprotocol && isunix,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_analysisprotocol);
  if forcecompute || todo,
    try
      fprintf(logfid,'AnalysisProtocol...\n');
      cmd = sprintf('./analysis_protocol.pl %s %s',expdir,analysis_protocol);
      unix(cmd);
    catch ME,
      msgs = {sprintf('Error running AnalysisProtocol:\n%s',getReport(ME))};
      fprintf(logfid,'AnalysisProtocol failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end; return;
    end
  end
  
  % make sure analysis protocol files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_analysisprotocol);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing analysisprotocol file %s',x),missingfiles,'UniformOutput',false);
      fprintf(logfid,'AnalysisProtocol failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    if logfid > 1, fclose(logfid); end; return;
  end
  
end


%% complete checks

stage = 'automaticchecks_complete';

if doautomaticcheckscomplete,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticcheckscomplete);
  if forcecompute || todo,
    try
      fprintf(logfid,'AutomaticChecks_Complete...\n');
      [success1,msgs] = FlyBowlAutomaticChecks_Complete(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        automaticcheckscomplete_params{:});
      
      if ~success1,
        fprintf(logfid,'AutomaticChecks_Complete failed:\n');
        fprintf(logfid,'%s\n',msgs{:});
        if logfid > 1, fclose(logfid); end; return;
      end
    catch ME,
      msgs = {sprintf('Error running AutomaticChecks_Complete:\n%s',getReport(ME))};
      fprintf(logfid,'AutomaticChecks_Complete failed:\n');
      fprintf(logfid,'%s\n',msgs{:});
      if logfid > 1, fclose(logfid); end; return;
    end
  end
  
  % make sure automatic checks files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticcheckscomplete);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing automaticchecks_complete file %s',x),missingfiles,'UniformOutput',false);
    fprintf(logfid,'AutomaticChecks_Complete failed:\n');
    fprintf(logfid,'%s\n',msgs{:});
    if logfid > 1, fclose(logfid); end; return;
  end
  
end

%% clean up

if isdeployed,
  delete(findall(0,'type','figure'));
end

fprintf(logfid,'Analysis pipeline completed!\n');
success = true;

if logfid > 1, fclose(logfid); end


function [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles)

ismissingfile = false;
missingfiles = {};
for i = 1:numel(requiredfiles),
  fn = requiredfiles{i};
  if isfield(dataloc_params,fn),
    fn = dataloc_params.(fn);
  end
  fn = fullfile(expdir,fn);
  if ~exist(fn,'file'),
    missingfiles{end+1} = fn; %#ok<AGROW>
    ismissingfile = true;
  end
end