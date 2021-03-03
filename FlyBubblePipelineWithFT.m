function [success,msgs,stage] = FlyBubblePipelineWithFT(expdir,varargin)
%runs the fly bubble pipeline 12/19/2016 - set up for Andy Seed's rig to
%run analysis on windows 

success = false;
msgs = {};
stage = 'start';

[analysis_protocol,...
 settingsdir,...
 datalocparamsfilestr,...
 automaticchecksincoming_params,...
 registration_params,...
 sexclassification_params,...
 computeperframefeatures_params,...
 computehoghofperframefeatures_params,...
 computeperframestats_params,...
 plotperframestats_params,...
 makectraxresultsmovie_params,...
 extradiagnostics_params,...
 automaticcheckscomplete_params,...
 forcecompute,...
 doautomaticchecksincoming,...
 doflytracking,...
 doregistration,...
 doledonoffdetection,...
 dosexclassification,...
 dotrackwings,...
 docomputeperframefeatures,...
 docomputehoghofperframefeatures,...
 dojaabadetect,...
 docomputeperframestats,...
 doplotperframestats,...
 domakectraxresultsmovie,...
 doextradiagnostics,...
 doanalysisprotocol,...
 doautomaticcheckscomplete,...
 requiredfiles_start,...
 requiredfiles_automaticchecksincoming,...
 requiredfiles_flytracker,...
 requiredfiles_registration,...
 requiredfiles_ledonoffdetection,...
 requiredfiles_sexclassification,...
 requiredfiles_wingtracking,...
 requiredfiles_computeperframefeatures,...
 requiredfiles_computehoghofperframefeatures,...
 requiredfiles_computeperframestats,...
 requiredfiles_plotperframestats,...
 requiredfiles_makectraxresultsmovie,...
 requiredfiles_extradiagnostics,...
 requiredfiles_analysisprotocol,...
 requiredfiles_automaticcheckscomplete] = ...
 myparse(varargin,...
         'analysis_protocol','current_bubble',...
         'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
         'datalocparamsfilestr','dataloc_params.txt',...
         'automaticchecksincoming_params',{},...
         'registration_params',{},...
         'sexclassification_params',{},...
         'computeperframefeatures_params',{},...
         'computeperframestats_params',{},...
         'computehoghofperframefeatures_params',{},...
         'plotperframestats_params',{},...
         'makectraxresultsmovie_params',{},...
         'extradiagnostics_params',{},...
         'automaticcheckscomplete_params',{},...
         'forcecompute',true,...
         'doautomaticchecksincoming',false,...
         'doflytracking',false, ...
         'doregistration',false,...
         'doledonoffdetection',false,...
         'dosexclassification',false,...
         'dotrackwings',false,...
         'docomputeperframefeatures',false,...
         'docomputehoghofperframefeatures',false,...
         'dojaabadetect',true,...
         'docomputeperframestats',false,...
         'doplotperframestats',false,...
         'domakectraxresultsmovie',true,...
         'doextradiagnostics',true,...
         'doanalysisprotocol',isunix,...
         'doautomaticcheckscomplete',true,...
         'requiredfiles_start',{},...
         'requiredfiles_automaticchecksincoming',{'automaticchecksincomingresultsfilestr'},...
         'requiredfiles_flytracker',{'ctraxfilestr'},...
         'requiredfiles_registration',{'trxfilestr','registrationmatfilestr','registrationtxtfilestr','registrationimagefilestr'},...
         'requiredfiles_ledonoffdetection',{'indicatordatafilestr'},...
         'requiredfiles_sexclassification',{'sexclassifierdiagnosticsfilestr'},...
         'requiredfiles_wingtracking',{'wingtrxfilestr'},...
         'requiredfiles_computeperframefeatures',{'perframedir', 'PERFRAMEMATFILES'},...         
         'requiredfiles_computehoghofperframefeatures',{'hs_sup_01_01_1.mat','hf_01_01_1.mat'},...
         'requiredfiles_computeperframestats',{'statsperframetxtfilestr','statsperframematfilestr','histperframetxtfilestr','histperframematfilestr'},...
         'requiredfiles_plotperframestats',{'figdir','ANALYSISPLOTS'},...
         'requiredfiles_makectraxresultsmovie',{'CTRAXRESULTSMOVIE'},...
         'requiredfiles_extradiagnostics',{'biasdiagnosticsimagefilestr','biasdiagnosticsfilestr','biasdiagnosticsmatfilestr',...
         'temperaturediagnosticsfilestr','bkgddiagnosticsimagefilestr','bkgddiagnosticsfilestr','bkgddiagnosticsmatfilestr',...
         'videodiagnosticsfilestr','videodiagnosticsmatfilestr','videodiagnosticsimagefilestr'},...
         'requiredfiles_analysisprotocol',{'analysis_protocol.txt'},...
         'requiredfiles_automaticcheckscomplete',{'automaticcheckscompleteresultsfilestr'});
%         'requiredfiles_start',{'annfilestr','ctraxfilestr'},...

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

if ischar(doledonoffdetection),
  doledonoffdetection = str2double(doledonoffdetection) ~= 0;
end
if ischar(dotrackwings),
  dotrackwings = str2double(dotrackwings) ~= 0;
end

if ischar(dosexclassification),
  dosexclassification = str2double(dosexclassification) ~= 0;
end

if ischar(docomputeperframefeatures),
  docomputeperframefeatures = str2double(docomputeperframefeatures) ~= 0;
end

if ischar(docomputehoghofperframefeatures),
  docomputehoghofperframefeatures = str2double(docomputehoghofperframefeatures) ~= 0;
end

if ischar(dojaabadetect),
  dojaabadetect = str2double(dojaabadetect) ~= 0;
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

% [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_start);
% if ismissingfile,
%   msgs = cellfun(@(x) sprintf('Missing start file %s',x),missingfiles,'UniformOutput',false);
%   fprintf('%s\n',msgs{:});
%   return;
% end



%% incoming checks
stage = 'automaticchecks_incoming';

if doautomaticchecksincoming,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticchecksincoming);
  if forcecompute || todo,
    fprintf('AutomaticChecks_Incoming...\n');

    [success1,msgs] = FlyBubbleAutomaticChecks_Incoming(expdir,...
      'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
      automaticchecksincoming_params{:});

    if ~success1,
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



%% Run FlyTracker
stage = 'flytracker' ;
if doflytracking ,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_flytracker) ;
  if forcecompute || todo ,    
    fprintf('FlyTracker...\n');      
    FlyTrackerWrapperForFlyBubble(expdir, settingsdir, analysis_protocol, dataloc_params, forcecompute) ;     
  end  
  
  % make sure flytracker files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_flytracker);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing FlyTracker file %s',x),missingfiles,'UniformOutput',false);
    fprintf('FlyTracker failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end



%% registration
stage = 'registration';
if doregistration,  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_registration);
  if forcecompute || todo,
    fprintf('RegisterTrx...\n');

    FlyBubbleRegisterTrx(expdir,...
      'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
      registration_params{:});
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



%% led indicator detection on/off
stage = 'ledonoffdetection';
if doledonoffdetection,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_ledonoffdetection);
  if forcecompute || todo, 
    fprintf('Detect LED onoff ...\n');
    FlyBubbleDectectIndicatorLedOnOff(expdir,...
      'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
      registration_params{:});
  end
  % make sure leddetection files exist requiredfiles_ledonoffdetection
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_ledonoffdetection);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing registration file %s',x),missingfiles,'UniformOutput',false);
    fprintf('RegisterTrx failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
end



%% Wing tracking and choose orientations
stage = 'trackwings';
if dotrackwings,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_wingtracking);
  if forcecompute || todo,
    fprintf('Wing tracking...\n');
    FlyTracker2WingTracking(expdir, ...
                            'dataloc_params', dataloc_params, ...
                            'settingsdir', settingsdir, ...
                            'analysis_protocol', analysis_protocol) ;
  end
  % make sure sexclassification files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_wingtracking);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing wingtracking file %s',x),missingfiles,'UniformOutput',false);
    fprintf('TrackWings failed:\n');
    fprintf('%s\n',msgs{:});
    return
  end
end

%% sex classification 

stage = 'sexclassification';

if dosexclassification,
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_sexclassification);
  if forcecompute || todo,    
    fprintf('ClassifySex...\n');
    FlyBubbleClassifySex(expdir,...
      'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
      sexclassification_params{:});
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
  
%% compute locomotor and social per-frame features

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
    fprintf('ComputePerFrameFeatures...\n');
    FlyBubbleComputePerFrameFeatures(expdir,...
      'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
      'forcecompute',forcecompute,...
      computeperframefeatures_params{:});
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

%% computer hog hof per-frame features
%flies_hoghof_hs_notaligned -- flow computed using Horn-Schunck. The current frame and next frame are not aligned in any way.

stage = 'computehoghofperframefeatures';

if docomputehoghofperframefeatures,
  requiredfiles_computehoghofperframefeatures = fullfile(dataloc_params.perframedir,requiredfiles_computehoghofperframefeatures);
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computehoghofperframefeatures);
  if forcecompute || todo,
    pwdprev = pwd;      
    datalocparamsfilestr = 'dataloc_params.txt';
    datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
    dataloc_params = ReadParams(datalocparamsfile);
    trxfilestr = fullfile(expdir,dataloc_params.trxfilestr);
    movfilestr = fullfile(expdir,dataloc_params.moviefilestr);
    spacetimefeaturesdir = fileparts(which('preparePerFrameFtrs'));
    cd (spacetimefeaturesdir);
    fprintf('preparePreFrameFtrs hog_hof...\n');
    preparePerFrameFtrs(movfilestr,trxfilestr,false,false);      
    cd(pwdprev);
  end
  % make sure computeperframefeatures files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframefeatures);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing computeperframefeatures file %s',x),missingfiles,'UniformOutput',false);
    fprintf('preparePreFrameFtrs hog_hof failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end

%% behavior detection
stage = 'jaabadetect';
if dojaabadetect,
    % For reasons I don't understand, this code errors if you remove the
    % try/catch... (ALT, 2020-02-27)
    try
        datalocparamsfilestr = 'dataloc_params.txt';
        datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
        dataloc_params = ReadParams(datalocparamsfile);
        jaabaclassifierparamsfilestrs = fullfile(settingsdir,analysis_protocol,dataloc_params.jaabaclassifierparamsfilestrs);
        jabfiles = textread(jaabaclassifierparamsfilestrs,'%s');
        fprintf('JAABAdetect...\n');
        JAABADetect(expdir, 'jabfiles', jabfiles, 'forcecompute', forcecompute, 'doteardownpath', true) ;
    catch ME,
        msgs = {sprintf('Error running JAABADetect:\n%s',getReport(ME))};
        fprintf('JAABADetect failed:\n');
        fprintf('%s\n',msgs{:});   
        return
    end    
end

    
%% make results movie

stage = 'ctraxresultsmovie';

if domakectraxresultsmovie,
  
  i = find(strcmp('CTRAXRESULTSMOVIE',requiredfiles_makectraxresultsmovie),1);
  if ~isempty(i),
    [~,basename] = fileparts(expdir);
    avifilestr = sprintf('%s_%s',dataloc_params.ctraxresultsavifilestr,basename);
    h264file = [avifilestr,'.mp4'];
    requiredfiles_makectraxresultsmovie{i} = h264file;
  end
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makectraxresultsmovie);
  if forcecompute || todo ,    
    fprintf('MakeCtraxResultsMovie...\n');
    FlyBubbleMakeCtraxResultsMovie(expdir,...
                                   'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
                                   makectraxresultsmovie_params{:});
  end
  
  % make sure makectraxresultsmovie files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makectraxresultsmovie);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing makectraxresultsmovie file %s',x),missingfiles,'UniformOutput',false);
    fprintf('MakeCtraxResultsMovie failed:\n');
    fprintf('%s\n',msgs{:});
    return
  end
  
end

%% complete checks

stage = 'automaticchecks_complete';

if doautomaticcheckscomplete,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticcheckscomplete);
  if forcecompute || todo,
    fprintf('AutomaticChecks_Complete...\n');
    [success1, msgs] = ...
      FlyBubbleAutomaticChecks_Complete(expdir,...
                                        'settingsdir',settingsdir, ...
                                        'analysis_protocol',analysis_protocol,...
                                        automaticcheckscomplete_params{:});
    if ~success1,
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

fprintf('Analysis pipeline completed!\n');
success = true;

