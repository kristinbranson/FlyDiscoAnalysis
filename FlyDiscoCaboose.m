function [success,msgs,stage] = FlyDiscoCaboose(expdir, varargin) 
% Runs just the auto-checks-complete from the FlyDisco pipeline.
  
% % Get info about the state of the repo, output to stdout
% this_script_path = mfilename('fullpath') ;
% source_folder_path = fileparts(this_script_path) ;
% git_report = get_git_report(source_folder_path) ;
% fprintf('%s', git_report) ;

% Initialize the outputs
success = false;
msgs = {};
stage = 'start';

% Process varargin parameters
if length(varargin)==1 && isstruct(varargin{1}) ,
  argument_parameters = varargin{1} ;
else
  argument_parameters = struct_from_name_value_list(varargin) ;
end

% First, get the settingsdir
settingsdir = lookup_in_struct(argument_parameters, 'settingsdir', '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings') ;        

% Get the metadata file path
metadata_file_path = determine_metadata_file_path(expdir) ;

% See if analysis protocol was passed in
try 
    analysis_protocol = argument_parameters.analysis_protocol ;
catch me
    if isequal(me.identifier, 'MATLAB:nonExistentField') ,
        % Read the experiment metadata to determine the analysis_protoocol
        % Also depends on what options exist in the settingsdir
        % Don't want to do this if analysis-protocol was specified in varargin,
        % because sometimes metadata doesn't specify analysis_protocol at all.
        analysis_protocol = analysis_protocol_from_metadata_file(metadata_file_path, settingsdir) ;
    else
        rethrow(me) ;
    end
end    

% Read in the analysis protocol parameters from the analysis-protocol folder, if
% it exists
analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
analysis_protocol_parameters_file_path = fullfile(analysis_protocol_folder_path, 'analysis-protocol-parameters.txt') ;
if exist(analysis_protocol_parameters_file_path, 'file') ,
    analysis_protocol_parameters = ReadParams(analysis_protocol_parameters_file_path) ;  % a struct
else
    analysis_protocol_parameters = struct() ;
end

% Set the default analysis parameters
default_analysis_parameters_as_list = ...
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
     'doautomaticchecksincoming','on',...
     'doflytracking','on', ...
     'doregistration','on',...
     'doledonoffdetection','on',...
     'dosexclassification','on',...
     'dotrackwings','off',...
     'docomputeperframefeatures','on',...
     'docomputehoghofperframefeatures','off',...
     'dojaabadetect','on',...
     'docomputeperframestats','on',...
     'doplotperframestats','on',...
     'domakectraxresultsmovie','on',...
     'doextradiagnostics','off',...
     'doanalysisprotocol','off',...
     'doautomaticcheckscomplete','off',...
     'requiredfiles_start',{},...
     'requiredfiles_automaticchecksincoming',{'automaticchecksincomingresultsfilestr'},...
     'requiredfiles_flytracker',{'ctraxfilestr'},...
     'requiredfiles_registration',{'trxfilestr','registrationmatfilestr','registrationtxtfilestr','registrationimagefilestr'},...
     'requiredfiles_ledonoffdetection',{'indicatordatafilestr'},...
     'requiredfiles_sexclassification',{'sexclassifierdiagnosticsfilestr'},...
     'requiredfiles_wingtracking',{'wingtrxfilestr'},...
     'requiredfiles_computeperframefeatures',{'perframedir', 'PERFRAMEMATFILES'},...
     'requiredfiles_computehoghofperframefeatures',{'hs_sup_01_01_1.mat','hf_01_01_1.mat'},...
     'requiredfiles_computeperframestats',{'statsperframetxtfilestr','statsperframematfilestr'},...
     'requiredfiles_plotperframestats',{'figdir','ANALYSISPLOTS'},...
     'requiredfiles_makectraxresultsmovie',{'CTRAXRESULTSMOVIE'},...
     'requiredfiles_extradiagnostics',{'biasdiagnosticsimagefilestr','biasdiagnosticsfilestr','biasdiagnosticsmatfilestr',...
     'temperaturediagnosticsfilestr','bkgddiagnosticsimagefilestr','bkgddiagnosticsfilestr','bkgddiagnosticsmatfilestr',...
     'videodiagnosticsfilestr','videodiagnosticsmatfilestr','videodiagnosticsimagefilestr'},...
     'requiredfiles_analysisprotocol',{'analysis_protocol.txt'},...
     'requiredfiles_automaticcheckscomplete',{'automaticcheckscompleteresultsfilestr'} } ;
 
 %put back in when hist is ready
 %      'requiredfiles_computeperframestats',{'statsperframetxtfilestr','statsperframematfilestr','histperframetxtfilestr','histperframematfilestr'},...
default_analysis_parameters = struct_from_name_value_list(default_analysis_parameters_as_list) ;
   
% Combine the default parameters with those from the analysis-protocol folder, giving precedence to the analysis-protocol folder ones 
analysis_parameters = merge_structs(default_analysis_parameters, analysis_protocol_parameters, argument_parameters) ;

% Assign the paramters to individual variables
datalocparamsfilestr = lookup_in_struct(analysis_parameters, 'datalocparamsfilestr') ;
automaticchecksincoming_params = lookup_in_struct(analysis_parameters, 'automaticchecksincoming_params') ;
registration_params = lookup_in_struct(analysis_parameters, 'registration_params') ;
sexclassification_params = lookup_in_struct(analysis_parameters, 'sexclassification_params') ;
computeperframefeatures_params = lookup_in_struct(analysis_parameters, 'computeperframefeatures_params') ;
makectraxresultsmovie_params = lookup_in_struct(analysis_parameters, 'makectraxresultsmovie_params') ;
automaticcheckscomplete_params = lookup_in_struct(analysis_parameters, 'automaticcheckscomplete_params') ;
doautomaticchecksincoming = lookup_in_struct(analysis_parameters, 'doautomaticchecksincoming') ;
doflytracking = lookup_in_struct(analysis_parameters, 'doflytracking') ;
doregistration = lookup_in_struct(analysis_parameters, 'doregistration') ;
doledonoffdetection = lookup_in_struct(analysis_parameters, 'doledonoffdetection') ;
dosexclassification = lookup_in_struct(analysis_parameters, 'dosexclassification') ;
dotrackwings = lookup_in_struct(analysis_parameters, 'dotrackwings') ;
docomputeperframefeatures = lookup_in_struct(analysis_parameters, 'docomputeperframefeatures') ;
docomputehoghofperframefeatures = lookup_in_struct(analysis_parameters, 'docomputehoghofperframefeatures') ;
dojaabadetect = lookup_in_struct(analysis_parameters, 'dojaabadetect') ;
docomputeperframestats = lookup_in_struct(analysis_parameters, 'docomputeperframestats') ;
doplotperframestats = lookup_in_struct(analysis_parameters, 'doplotperframestats') ;
domakectraxresultsmovie = lookup_in_struct(analysis_parameters, 'domakectraxresultsmovie') ;
doextradiagnostics = lookup_in_struct(analysis_parameters, 'doextradiagnostics') ;
doanalysisprotocol = lookup_in_struct(analysis_parameters, 'doanalysisprotocol') ; 
doautomaticcheckscomplete = lookup_in_struct(analysis_parameters, 'doautomaticcheckscomplete') ;
requiredfiles_automaticchecksincoming = lookup_in_struct(analysis_parameters, 'requiredfiles_automaticchecksincoming') ;
requiredfiles_flytracker = lookup_in_struct(analysis_parameters, 'requiredfiles_flytracker') ;
requiredfiles_registration = lookup_in_struct(analysis_parameters, 'requiredfiles_registration') ;
requiredfiles_ledonoffdetection = lookup_in_struct(analysis_parameters, 'requiredfiles_ledonoffdetection') ;
requiredfiles_sexclassification = lookup_in_struct(analysis_parameters, 'requiredfiles_sexclassification') ;
requiredfiles_wingtracking = lookup_in_struct(analysis_parameters, 'requiredfiles_wingtracking') ;
requiredfiles_computeperframefeatures = lookup_in_struct(analysis_parameters, 'requiredfiles_computeperframefeatures') ;
requiredfiles_computeperframestats = lookup_in_struct(analysis_parameters, 'requiredfiles_computeperframestats') ;
requiredfiles_computehoghofperframefeatures = lookup_in_struct(analysis_parameters, 'requiredfiles_computehoghofperframefeatures') ;
requiredfiles_makectraxresultsmovie = lookup_in_struct(analysis_parameters, 'requiredfiles_makectraxresultsmovie') ;
requiredfiles_automaticcheckscomplete = lookup_in_struct(analysis_parameters, 'requiredfiles_automaticcheckscomplete') ;

% Read in the dataloc params         
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

% Coerse all the do* variables to on/off/force
doautomaticchecksincoming = coerce_to_on_off_force(doautomaticchecksincoming) ;
doflytracking = coerce_to_on_off_force(doflytracking) ;
doregistration = coerce_to_on_off_force(doregistration) ;
doledonoffdetection = coerce_to_on_off_force(doledonoffdetection) ;
dosexclassification = coerce_to_on_off_force(dosexclassification) ;
dotrackwings = coerce_to_on_off_force(dotrackwings) ;
docomputeperframefeatures = coerce_to_on_off_force(docomputeperframefeatures) ;
docomputehoghofperframefeatures = coerce_to_on_off_force(docomputehoghofperframefeatures) ;
dojaabadetect = coerce_to_on_off_force(dojaabadetect) ;
docomputeperframestats = coerce_to_on_off_force(docomputeperframestats) ; %#ok<NASGU>
doplotperframestats = coerce_to_on_off_force(doplotperframestats) ; %#ok<NASGU>
domakectraxresultsmovie = coerce_to_on_off_force(domakectraxresultsmovie) ;
doextradiagnostics = coerce_to_on_off_force(doextradiagnostics) ; %#ok<NASGU>
doanalysisprotocol = coerce_to_on_off_force(doanalysisprotocol) ; %#ok<NASGU>
doautomaticcheckscomplete = coerce_to_on_off_force(doautomaticcheckscomplete) ;

% Print the settings values in use
variable_names_to_print = ...
    { 'settingsdir', ...
      'analysis_protocol', ...
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

% % Print the canonical path to the analysis folder
% canonical_analysis_protocol_folder_path = realpath(analysis_protocol_folder_path) ;
% fprintf('Canonical path to analysis protocol folder is:\n  %s\n\n', canonical_analysis_protocol_folder_path) ;

%% check that experiment exists
if ~exist(expdir,'dir'),
  msgs = {sprintf('Experiment directory %s does not exist',expdir)};
  fprintf('%s\n',msgs{:});
  return;
end


%% automaticchecks_complete
% It's the only thing we do.
stage = 'automaticchecks_complete';
if is_on_or_force(doautomaticcheckscomplete) ,
  forcecompute = is_force(doautomaticcheckscomplete) ;  
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

success = true;

