function [success,msgs,stage] = FlyBowlAnalysisPipeline(expdir,varargin)

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
  docomputeperframefeatures,docomputeperframestats,doplotperframestats,...
  domakectraxresultsmovie,doextradiagnostics,doanalysisprotocol,...
  doautomaticcheckscomplete,...
  requiredfiles_start,requiredfiles_automaticchecksincoming,...
  requiredfiles_registration,...
  requiredfiles_sexclassification,...
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
  'docomputeperframefeatures',true,...
  'docomputeperframestats',true,...
  'doplotperframestats',true,...
  'domakectraxresultsmovie',true,...
  'doextradiagnostics',true,...
  'doanalysisprotocol',isunix,...
  'doautomaticcheckscomplete',true,...
  'requiredfiles_start',{'annfilestr','ctraxfilestr'},...
  'requiredfiles_automaticchecksincoming',{'automaticchecksincomingresultsfilestr'},...
  'requiredfiles_registration',{'trxfilestr','registrationmatfilestr','registrationtxtfilestr','registrationimagefilestr'},...
  'requiredfiles_sexclassification',{'sexclassifierdiagnosticsfilestr'},...
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

if ischar(forcecompute),
  forcecompute = str2double(forcecompute) ~= 0;
end

if ischar(doautomaticchecksincoming),
  doautomaticchecksincoming = str2double(doautomaticchecksincoming) ~= 0;
end

if ischar(doregistration),
  doregistration = str2double(doregistration) ~= 0;
end

if ischar(dosexclassification),
  dosexclassification = str2double(dosexclassification) ~= 0;
end

if ischar(docomputeperframefeatures),
  docomputeperframefeatures = str2double(docomputeperframefeatures) ~= 0;
end

if ischar(docomputeperframestats),
  docomputeperframestats = str2double(docomputeperframestats) ~= 0;
end

if ischar(doplotperframestats),
  doplotperframestats = str2double(doplotperframestats) ~= 0;
end

if ischar(domakectraxresultsmovie),
  domakectraxresultsmovie = str2double(domakectraxresultsmovie) ~= 0;
end

if ischar(doextradiagnostics),
  doextradiagnostics = str2double(doextradiagnostics) ~= 0;
end

if ischar(doautomaticcheckscomplete),
  doautomaticcheckscomplete = str2double(doautomaticcheckscomplete) ~= 0;
end

%% check that experiment exists

if ~exist(expdir,'dir'),
  msgs = {sprintf('Experiment directory %s does not exist',expdir)};
  fprintf('%s\n',msgs{:});
  return;
end

%% check for files required at start

[ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_start);
if ismissingfile,
  msgs = cellfun(@(x) sprintf('Missing start file %s',x),missingfiles,'UniformOutput',false);
  fprintf('%s\n',msgs{:});
  return;
end

%% incoming checks

stage = 'automaticchecks_incoming';

if doautomaticchecksincoming,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticchecksincoming);
  if forcecompute || todo,
    try
      fprintf('AutomaticChecks_Incoming...\n');

      [success1,msgs] = FlyBowlAutomaticChecks_Incoming(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        automaticchecksincoming_params{:});
      
      if ~success1,
        fprintf('AutomaticChecks_Incoming failed:\n');
        fprintf('%s\n',msgs{:});
        return;
      end
    catch ME,
      msgs = {sprintf('Error running AutomaticChecks_Incoming:\n%s',getReport(ME))};
      fprintf('AutomaticChecks_Incoming failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure automatic checks files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticchecksincoming);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing automaticchecks_incoming file %s',x),missingfiles,'UniformOutput',false);
    fprintf('AutomaticChecks_Incoming failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end

%% registration

stage = 'registration';

if doregistration,
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_registration);
  if forcecompute || todo,
    
    try
      fprintf('RegisterTrx...\n');

      FlyBowlRegisterTrx(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        registration_params{:});
      
    catch ME,
      msgs = {sprintf('Error running RegisterTrx:\n%s',getReport(ME))};
      fprintf('RegisterTrx failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure registration files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_registration);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing registration file %s',x),missingfiles,'UniformOutput',false);
    fprintf('RegisterTrx failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end

%% sex classification


stage = 'sexclassification';

if dosexclassification,
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_sexclassification);
  if forcecompute || todo,
    
    try
      fprintf('ClassifySex...\n');
      FlyBowlClassifySex2(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        sexclassification_params{:});
      
    catch ME,
      msgs = {sprintf('Error running SexClassification:\n%s',getReport(ME))};
      fprintf('SexClassification failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure sexclassification files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_sexclassification);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing sexclassification file %s',x),missingfiles,'UniformOutput',false);
    fprintf('SexClassification failed:\n');
    fprintf('%s\n',msgs{:});
    return;
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
      fprintf('ComputePerFrameFeatures...\n');
      FlyBowlComputePerFrameFeatures(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        'forcecompute',forcecompute,...
        computeperframefeatures_params{:});
      
    catch ME,
      msgs = {sprintf('Error running ComputePerFrameFeatures:\n%s',getReport(ME))};
      fprintf('ComputePerFrameFeatures failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure computeperframefeatures files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframefeatures);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing computeperframefeatures file %s',x),missingfiles,'UniformOutput',false);
    fprintf('ComputePerFrameFeatures failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end

%% compute per-frame statistics

stage = 'computeperframestats';

if docomputeperframestats,
    
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframestats);
  if forcecompute || todo,
    
    try
      fprintf('ComputePerFrameStats...\n');
      FlyBowlComputePerFrameStats2(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        computeperframestats_params{:});
      
    catch ME,
      msgs = {sprintf('Error running ComputePerFrameStats:\n%s',getReport(ME))};
      fprintf('ComputePerFrameStats failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure computeperframestats files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframestats);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing computeperframestats file %s',x),missingfiles,'UniformOutput',false);
    fprintf('ComputePerFrameStats failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end

%% plot per-frame stats

stage = 'plotperframestats';

if doplotperframestats,
    
  if ismember('ANALYSISPLOTS',requiredfiles_plotperframestats),
    histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
    hist_perframefeatures = ReadHistPerFrameFeatures(histperframefeaturesfile);
    hist_fields = unique({hist_perframefeatures.field});
    for i = 1:numel(hist_fields),
      savename = sprintf('hist_%s.png',hist_fields{i});
      savename = fullfile(dataloc_params.figdir,savename);
      requiredfiles_plotperframestats{end+1} = savename; %#ok<AGROW>
    end
    requiredfiles_plotperframestats = setdiff(requiredfiles_plotperframestats,{'ANALYSISPLOTS'});
    requiredfiles_plotperframestats{end+1} = fullfile(dataloc_params.figdir,'stats.png');
  end
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_plotperframestats);
  if forcecompute || todo,
    
    try
      fprintf('PlotPerFrameStats...\n');
      FlyBowlPlotPerFrameStats2(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        plotperframestats_params{:});
    catch ME,
      msgs = {sprintf('Error running PlotPerFrameStats:\n%s',getReport(ME))};
      fprintf('PlotPerFrameStats failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure plotperframestats files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_plotperframestats);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing plotperframestats file %s',x),missingfiles,'UniformOutput',false);
    fprintf('PlotPerFrameStats failed:\n');
    fprintf('%s\n',msgs{:});
    return;
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
      fprintf('MakeCtraxResultsMovie...\n');
      FlyBowlMakeCtraxResultsMovie(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        makectraxresultsmovie_params{:});
    catch ME,
      msgs = {sprintf('Error running MakeCtraxResultsMovie:\n%s',getReport(ME))};
      fprintf('MakeCtraxResultsMovie failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure makectraxresultsmovie files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makectraxresultsmovie);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing makectraxresultsmovie file %s',x),missingfiles,'UniformOutput',false);
    fprintf('MakeCtraxResultsMovie failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end

%% extra diagnostics

stage = 'extradiagnostics';

if doextradiagnostics,
    
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_extradiagnostics);
  if forcecompute || todo,
    
    try
      fprintf('ExtraDiagnostics...\n');
      FlyBowlExtraDiagnostics(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        extradiagnostics_params{:});
    catch ME,
      msgs = {sprintf('Error running ExtraDiagnostics:\n%s',getReport(ME))};
      fprintf('ExtraDiagnostics failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure extradiagnostics files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_extradiagnostics);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing extradiagnostics file %s',x),missingfiles,'UniformOutput',false);
    fprintf('ExtraDiagnostics failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end

%% generate analysis protocol file

stage = 'analysisprotocol';

if doanalysisprotocol && isunix,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_analysisprotocol);
  if forcecompute || todo,
    try
      fprintf('AnalysisProtocol...\n');
      cmd = sprintf('./analysis_protocol.pl %s %s',expdir,analysis_protocol);
      unix(cmd);
    catch ME,
      msgs = {sprintf('Error running AnalysisProtocol:\n%s',getReport(ME))};
      fprintf('AnalysisProtocol failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure analysis protocol files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_analysisprotocol);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing analysisprotocol file %s',x),missingfiles,'UniformOutput',false);
      fprintf('AnalysisProtocol failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end


%% complete checks

stage = 'automaticchecks_complete';

if doautomaticcheckscomplete,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticcheckscomplete);
  if forcecompute || todo,
    try
      fprintf('AutomaticChecks_Complete...\n');
      [success1,msgs] = FlyBowlAutomaticChecks_Complete(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        automaticcheckscomplete_params{:});
      
      if ~success1,
        fprintf('AutomaticChecks_Complete failed:\n');
        fprintf('%s\n',msgs{:});
        return;
      end
    catch ME,
      msgs = {sprintf('Error running AutomaticChecks_Complete:\n%s',getReport(ME))};
      fprintf('AutomaticChecks_Complete failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure automatic checks files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticcheckscomplete);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing automaticchecks_complete file %s',x),missingfiles,'UniformOutput',false);
    fprintf('AutomaticChecks_Complete failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end

%%

fprintf('Analysis pipeline completed!\n');
success = true;

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