function [success,msgs,stage] = FlyDiscoPipeline(expdir, varargin)
%runs the fly bubble pipeline 12/19/2016 - set up for Andy Seed's rig to
%run analysis on windows 

% Get info about the state of the repo, output to stdout
this_script_path = mfilename('fullpath') ;
source_folder_path = fileparts(this_script_path) ;
git_report = get_git_report(source_folder_path) ;
fprintf('%s', git_report) ;

% Initialize the outputs
success = false;
msgs = {};
stage = 'start';

% First, get the settingsdir
settingsdir = lookup_in_name_value_list(varargin, 'settingsdir', '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings') ;
        
% Read the experiment metadata to determine the analysis_protoocol
% Also depends on what options exist in the settingsdir
metadata_file_path = determine_metadata_file_path(expdir) ;
analysis_protocol_from_experiment_metadata = analysis_protocol_from_metadata_file(metadata_file_path, settingsdir) ;    

% Allow the varargins to override the analysis_protocol
analysis_protocol = ...
    lookup_in_name_value_list(varargin,...
                              'analysis_protocol', analysis_protocol_from_experiment_metadata) ;
        
% Read in the analysis protocol parameters from the analysis-protocol folder, if
% it exists
analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
analysis_protocol_parameters_file_path = fullfile(analysis_protocol_folder_path, 'analysis-protocol-parameters.yaml') ;
if exist(analysis_protocol_parameters_file_path, 'file') ,
    analysis_protocol_parameters = parse_simple_yaml_file(analysis_protocol_parameters_file_path) ;
else
    analysis_protocol_parameters = cell(1,0) ;
end

% Set the default analysis parameters
default_analysis_parameters = ...
    {'datalocparamsfilestr','dataloc_params.txt',...
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
     'doautomaticchecksincoming',true,...
     'doflytracking',true, ...
     'doregistration',true,...
     'doledonoffdetection',true,...
     'dosexclassification',true,...
     'dotrackwings',false,...
     'docomputeperframefeatures',true,...
     'docomputehoghofperframefeatures',false,...
     'dojaabadetect',true,...
     'docomputeperframestats',true,...
     'doplotperframestats',true,...
     'domakectraxresultsmovie',true,...
     'doextradiagnostics',true,...
     'doanalysisprotocol',isunix(),...
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
     'requiredfiles_automaticcheckscomplete',{'automaticcheckscompleteresultsfilestr'} } ;
     
% Combine the default parameters with those from the analysis-protocol folder, giving precedence to the analysis-protocol folder ones 
analysis_parameters = merge_name_value_lists(default_analysis_parameters, analysis_protocol_parameters) ;
     
% Assign the paramters to individual variables
datalocparamsfilestr = lookup_in_name_value_list(analysis_parameters, 'datalocparamsfilestr') ;
automaticchecksincoming_params = lookup_in_name_value_list(analysis_parameters, 'automaticchecksincoming_params') ;
registration_params = lookup_in_name_value_list(analysis_parameters, 'registration_params') ;
sexclassification_params = lookup_in_name_value_list(analysis_parameters, 'sexclassification_params') ;
computeperframefeatures_params = lookup_in_name_value_list(analysis_parameters, 'computeperframefeatures_params') ;
makectraxresultsmovie_params = lookup_in_name_value_list(analysis_parameters, 'makectraxresultsmovie_params') ;
automaticcheckscomplete_params = lookup_in_name_value_list(analysis_parameters, 'automaticcheckscomplete_params') ;
forcecompute = lookup_in_name_value_list(analysis_parameters, 'forcecompute') ;
doautomaticchecksincoming = lookup_in_name_value_list(analysis_parameters, 'doautomaticchecksincoming') ;
doflytracking = lookup_in_name_value_list(analysis_parameters, 'doflytracking') ;
doregistration = lookup_in_name_value_list(analysis_parameters, 'doregistration') ;
doledonoffdetection = lookup_in_name_value_list(analysis_parameters, 'doledonoffdetection') ;
dosexclassification = lookup_in_name_value_list(analysis_parameters, 'dosexclassification') ;
dotrackwings = lookup_in_name_value_list(analysis_parameters, 'dotrackwings') ;
docomputeperframefeatures = lookup_in_name_value_list(analysis_parameters, 'docomputeperframefeatures') ;
docomputehoghofperframefeatures = lookup_in_name_value_list(analysis_parameters, 'docomputehoghofperframefeatures') ;
dojaabadetect = lookup_in_name_value_list(analysis_parameters, 'dojaabadetect') ;
docomputeperframestats = lookup_in_name_value_list(analysis_parameters, 'docomputeperframestats') ;
doplotperframestats = lookup_in_name_value_list(analysis_parameters, 'doplotperframestats') ;
domakectraxresultsmovie = lookup_in_name_value_list(analysis_parameters, 'domakectraxresultsmovie') ;
doextradiagnostics = lookup_in_name_value_list(analysis_parameters, 'doextradiagnostics') ;
doanalysisprotocol = lookup_in_name_value_list(analysis_parameters, 'doanalysisprotocol') ; %#ok<NASGU>
doautomaticcheckscomplete = lookup_in_name_value_list(analysis_parameters, 'doautomaticcheckscomplete') ;
requiredfiles_automaticchecksincoming = lookup_in_name_value_list(analysis_parameters, 'requiredfiles_automaticchecksincoming') ;
requiredfiles_flytracker = lookup_in_name_value_list(analysis_parameters, 'requiredfiles_flytracker') ;
requiredfiles_registration = lookup_in_name_value_list(analysis_parameters, 'requiredfiles_registration') ;
requiredfiles_ledonoffdetection = lookup_in_name_value_list(analysis_parameters, 'requiredfiles_ledonoffdetection') ;
requiredfiles_sexclassification = lookup_in_name_value_list(analysis_parameters, 'requiredfiles_sexclassification') ;
requiredfiles_wingtracking = lookup_in_name_value_list(analysis_parameters, 'requiredfiles_wingtracking') ;
requiredfiles_computeperframefeatures = lookup_in_name_value_list(analysis_parameters, 'requiredfiles_computeperframefeatures') ;
requiredfiles_computehoghofperframefeatures = lookup_in_name_value_list(analysis_parameters, 'requiredfiles_computehoghofperframefeatures') ;
requiredfiles_makectraxresultsmovie = lookup_in_name_value_list(analysis_parameters, 'requiredfiles_makectraxresultsmovie') ;
requiredfiles_automaticcheckscomplete = lookup_in_name_value_list(analysis_parameters, 'requiredfiles_automaticcheckscomplete') ;

% Read in the dataloc params         
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
  docomputeperframestats = str2double(docomputeperframestats) ~= 0;  %#ok<NASGU>
end

if ischar(doplotperframestats),
  doplotperframestats = str2double(doplotperframestats) ~= 0;  %#ok<NASGU>
end

if ischar(domakectraxresultsmovie),
  domakectraxresultsmovie = str2double(domakectraxresultsmovie) ~= 0;
end

if ischar(doextradiagnostics),
  doextradiagnostics = str2double(doextradiagnostics) ~= 0;  %#ok<NASGU>
end

if ischar(doautomaticcheckscomplete),
  doautomaticcheckscomplete = str2double(doautomaticcheckscomplete) ~= 0;
end

% Print the settings values in use
variable_names_to_print = ...
    { 'settingsdir', ...
      'analysis_protocol', ...
      'forcecompute', ...
      'doautomaticchecksincoming', ...
      'doflytracking',  ...
      'doregistration', ...
      'doledonoffdetection', ...
      'dosexclassification', ...
      'dotrackwings', ...
      'docomputeperframefeatures', ...
      'docomputehoghofperframefeatures', ...
      'dojaabadetect', ...
      'docomputeperframestats', ...
      'doplotperframestats', ...
      'domakectraxresultsmovie', ...
      'doextradiagnostics', ...
      'doanalysisprotocol', ...
      'doautomaticcheckscomplete' }' ;
fprintf('Settings values in FlyDiscoPipeline():\n') ;    
for i = 1 : length(variable_names_to_print) ,
  variable_name = variable_names_to_print{i} ;
  value = eval(variable_name) ;
  fprintf('  %s: %s\n', variable_name, char_array_from_value(value)) ;
end 
fprintf('\n') ;

% Print the canonical path to the analysis folder
escaped_analysis_protocol_folder_path = escape_path_for_bash(analysis_protocol_folder_path) ;
command_line = sprintf('realpath %s', escaped_analysis_protocol_folder_path) ;
stdout = system_with_error_handling(command_line) ;
canonical_analysis_protocol_folder_path = strtrim(stdout) ;
fprintf('Canonical path to analysis protocol folder is:\n  %s\n\n', canonical_analysis_protocol_folder_path) ;

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
    fprintf('Running incoming automatic checks...\n');

    [success1,msgs] = FlyDiscoAutomaticChecksIncoming(expdir,...
      'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
      automaticchecksincoming_params{:});

    if ~success1,
      fprintf('Incoming automatic checks failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure automatic checks files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticchecksincoming);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing incoming automatic checks file %s',x),missingfiles,'UniformOutput',false);
    fprintf('Incoming automatic checks failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end



%% Run FlyTracker
stage = 'flytracker' ;
if doflytracking ,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_flytracker) ;
  if forcecompute || todo ,    
    fprintf('Running FlyTracker...\n');      
    FlyTrackerWrapperForFlyDisco(expdir, settingsdir, analysis_protocol, dataloc_params, forcecompute) ;     
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
    fprintf('Registering tracks...\n');

    FlyDiscoRegisterTrx(expdir,...
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
    fprintf('Detecting LED on/off transitions...\n');
    FlyDiscoDectectIndicatorLedOnOff(expdir,...
      'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
      registration_params{:});
  end
  % make sure leddetection files exist requiredfiles_ledonoffdetection
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_ledonoffdetection);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing registration file %s',x),missingfiles,'UniformOutput',false);
    fprintf('Detecting LED on/off transitions failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
end



%% Wing tracking and choose orientations
stage = 'trackwings';
if dotrackwings,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_wingtracking);
  if forcecompute || todo,
    fprintf('Running wing tracking...\n');
    FlyTracker2WingTracking(expdir, ...
                            'dataloc_params', dataloc_params, ...
                            'settingsdir', settingsdir, ...
                            'analysis_protocol', analysis_protocol) ;
  end
  % make sure sexclassification files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_wingtracking);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing wing tracking file %s',x),missingfiles,'UniformOutput',false);
    fprintf('Tracking wings failed:\n');
    fprintf('%s\n',msgs{:});
    return
  end
end

%% sex classification 

stage = 'sexclassification';

if dosexclassification,
  
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_sexclassification);
  if forcecompute || todo,    
    fprintf('Running sex classification...\n');
    FlyDiscoClassifySex(expdir,...
      'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
      sexclassification_params{:});
  end
  
  % make sure sexclassification files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_sexclassification);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing sex classification file %s',x),missingfiles,'UniformOutput',false);
    fprintf('Sex classification failed:\n');
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
    fprintf('Computing per-frame features...\n');
    FlyDiscoComputePerFrameFeatures(expdir,...
      'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
      'forcecompute',forcecompute,...
      computeperframefeatures_params{:});
  end
  
  % make sure computeperframefeatures files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframefeatures);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing per-frame features file %s',x),missingfiles,'UniformOutput',false);
    fprintf('Computing per-frame features failed:\n');
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
    fprintf('Computing HOG/HOF per-frame features...\n');
    preparePerFrameFtrs(movfilestr,trxfilestr,false,false);      
    cd(pwdprev);
  end
  % make sure computeperframefeatures files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_computeperframefeatures);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing HOG/HOF per-frame features file %s',x),missingfiles,'UniformOutput',false);
    fprintf('Computing HOG/HOF per-frame features failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end
  
end

%% behavior detection
%stage = 'jaabadetect'; 
if dojaabadetect,
  JAABADetectWrapper(expdir, settingsdir, analysis_protocol, forcecompute) ;
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
    fprintf('Making Ctrax results movie...\n');
    FlyDiscoMakeCtraxResultsMovie(expdir,...
                                   'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
                                   makectraxresultsmovie_params{:});
  end
  
  % make sure makectraxresultsmovie files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_makectraxresultsmovie);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing Ctrax results movie file %s',x),missingfiles,'UniformOutput',false);
    fprintf('Making Ctrax results movie failed:\n');
    fprintf('%s\n',msgs{:});
    return
  end
  
end

%% complete checks

stage = 'automaticchecks_complete';

if doautomaticcheckscomplete,
  todo = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticcheckscomplete);
  if forcecompute || todo,
    fprintf('Running completion automatic checks...\n');
    [success1, msgs] = ...
      FlyDiscoAutomaticChecksComplete(expdir,...
                                        'settingsdir',settingsdir, ...
                                        'analysis_protocol',analysis_protocol,...
                                        automaticcheckscomplete_params{:});
    if ~success1,
      fprintf('Running completion automatic checks failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure automatic checks files exist
  [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles_automaticcheckscomplete);
  if ismissingfile,
    msgs = cellfun(@(x) sprintf('Missing completion automatic checks file %s',x),missingfiles,'UniformOutput',false);
    fprintf('completion automatic checks failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end  
end

% If get here, analysis has completed successfully
fprintf('Analysis pipeline completed!\n');
success = true;

