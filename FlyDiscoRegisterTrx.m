function FlyDiscoRegisterTrx(expdir, varargin)
% Read in the trajectories in pixel units, save them out again augmented with
% trajectories in physical units.  Do this, have to first find the registration
% mark(s) in a 'typical' movie frame.

% Need this 'global' in several places.
trx_file_names_that_are_not_per_frame = ...
  {'id','moviename','annname','firstframe','arena','off',...
   'nframes','endframe','matname','fps','pxpermm'} ;

%
% Parse arguments
%
[analysis_protocol,settingsdir,registrationparamsfilestr,datalocparamsfilestr,dotemporalreg,dotemporaltruncation] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir',default_settings_folder_path(),...
  'registrationparamsfilestr','registration_params.txt',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'dotemporalreg',false,'dotemporaltruncation',false);


%
% Read in the data file locations
%
analysis_protocol_folder_path = fullfile(settingsdir,analysis_protocol) ;
datalocparamsfile = fullfile(analysis_protocol_folder_path,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);


%
% Read in the registration params
%

% name of parameters file
registrationparamsfile = fullfile(analysis_protocol_folder_path,registrationparamsfilestr);
if ~exist(registrationparamsfile,'file'),
  error('Registration params file %s does not exist',registrationparamsfile);
end

% read the file contents into a struct
registration_params_0 = ReadParams(registrationparamsfile);

% Determine final dotemporalreg, dotemporaltruncation
[dotemporalreg, dotemporaltruncation, idealVideoDuration] = determineTemporalStuff(registration_params_0, dotemporalreg, dotemporaltruncation) ;
if dotemporalreg && dotemporaltruncation ,
  error('Temporal registration and temporal truncation are mutually exclusive.  You can''t do both!') ;
end
if dotemporaltruncation && isempty(idealVideoDuration) ,
  error('Temporal truncation requested, but the idealVideoDuration is []') ;
end

%
% Detect the registration mark(s)
%

% Load the background file output by FlyTracker
if isfield(dataloc_params,'flytrackerbgstr')
  flytrackerbgfile = fullfile(expdir,dataloc_params.flytrackerbgstr);
  load(flytrackerbgfile,'bg');
  bg_mean = 255*bg.bg_mean;
else
  % Note: We really do want to error if this is missing.
  % We've been bitten by this not being what we thought it was.
  error('dataloc_params is missing field flytrackerbgstr, which is required');
end

% Load metadata
metadata = collect_metadata(expdir, dataloc_params.metadatafilestr) ;

% Tweak registration_params.bowlMarkerType, if present
if isfield(registration_params_0,'bowlMarkerType'),
  bowlMarkerType = determineBowlMarkerType(registration_params_0.bowlMarkerType, metadata, analysis_protocol_folder_path) ;  
  registration_params_1 = setfield(registration_params_0, 'bowlMarkerType', bowlMarkerType) ;  %#ok<SFLD> 
else
  registration_params_1 = registration_param_0 ;
end

% Tweak registration_params.maxDistCornerFrac_BowlLabel, if present
if isfield(registration_params_1,'maxDistCornerFrac_BowlLabel') ,
  maxDistCornerFrac_BowlLabel = ...
    determineMaxDistCornerFracLabel(registration_params_1.maxDistCornerFrac_BowlLabel, metadata, 'maxDistCornerFrac_BowlLabel') ;
  registration_params_2 = setfield(registration_params_1, 'maxDistCornerFrac_BowlLabel', maxDistCornerFrac_BowlLabel) ;  %#ok<SFLD> 
else
  registration_params_2 = registration_param_1 ;
end

% Call the core registration mark detection routine
registration_params_cell = marshallRegistrationParamsForDetectRegistrationMarks(registration_params_2, expdir, dataloc_params.registrationimagefilestr) ;
registration_data_0 = detectRegistrationMarks(registration_params_cell{:},'bkgdImage',bg_mean,'useNormXCorr',true);
fprintf('Detected registration marks.\n');


%
% Apply spatial registration to the trajectories
%

% Get name of input trx mat file
ctraxfile = fullfile(expdir,dataloc_params.ctraxfilestr);

% Get name of movie
moviefile = fullfile(expdir,dataloc_params.moviefilestr);

% Load trajectories
[trx,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile,'annname','');
if ~succeeded,
  error('Could not load trajectories from file %s',ctraxfile);
end

% Postprocess trajectories to remove nans from flytracker outputs
nids0 = numel(trx);
if ~isfield(registration_params_2,'maxFlyTrackerNanInterpFrames'),
  fprintf('maxFlyTrackerNanInterpFrames not set in registration_params, using default value.\n');
  args = {};
else
  args = {'maxFlyTrackerNanInterpFrames',registration_params_2.maxFlyTrackerNanInterpFrames};
end
[trx,ninterpframes,newid2oldid] = PostprocessFlyTracker(trx,args{:});
fprintf('Removed nans from tracker output.\n');
fprintf('Number of nans interpolated through: %d frames\n',ninterpframes);
fprintf('Number of identities was %d, now %d\n',nids0,numel(trx));

% Store some metadata about nan-removal in registration_data
registration_data_0p5 = setfield(registration_data_0, 'flytracker_nnanframes', ninterpframes) ;  %#ok<SFLD> 
registration_data_1 = setfield(registration_data_0p5, 'flytracker_nids0', nids0) ;  %#ok<SFLD> 

% Compute the median frame interval, if called for
if registration_params_2.usemediandt,
  meddt = medianIntervalFromTimestamps(timestamps) ;
else
  meddt = [] ;
end

% If there are zero tracks, something is wrong
if isempty(trx),
  error('No flies tracked.');
end

% Apply spatial registration to trajectories
trx = appendPhysicalUnitFieldsToTrx(trx, registration_data_1, meddt) ;
fprintf('Applied spatial registration.\n') ;


%
% Handle all the stuff that is specific to optogenetic experiments
%
registration_data_2 = ...
  handleOptogeneticExperimentRegistration(registration_data_1, registration_params_2, ...
                                          metadata, expdir, dataloc_params, timestamps, analysis_protocol_folder_path) ;
% registation_data_2 has the field ledIndicatorPoints, which
% registration_data_1 lacks.

%
% Crop start and end of trajectories based on fly loaded time
%
[trx, timestamps, registration_data_3, newid2oldid] = ...
  performTemporalRegistration(dotemporalreg, trx, timestamps, registration_data_2, newid2oldid, ...
                              trx_file_names_that_are_not_per_frame, expdir, dataloc_params) ;


% 
% Truncate end of movie based on value from registration params 
%
[trx, registration_data] = ...
  performTemporalTruncation(dotemporaltruncation, trx, registration_data_3, newid2oldid, ...
                            trx_file_names_that_are_not_per_frame, idealVideoDuration, moviefile) ;


%
% Save registered trx to file
%
trxfile = fullfile(expdir,dataloc_params.trxfilestr);
try
  if exist(trxfile,'file'),
    delete(trxfile);
  end
  save(trxfile,'trx','timestamps'); %??? which timestamps should this be? Shouldn't it be truncated to match registration data? 
  fprintf('Saved registered trx to file %s\n',trxfile);
catch exception
  fprintf('Could not save registered trx:\n%s\n',getReport(exception));
end


%
% Save registration_data to mat file
%
registrationmatfile = fullfile(expdir,dataloc_params.registrationmatfilestr);
registration_data_without_registerfn = rmfield(registration_data,'registerfn'); 
try
  if exist(registrationmatfile,'file'),
    delete(registrationmatfile);
  end
  save(registrationmatfile,'-struct','registration_data_without_registerfn');
  fprintf('Saved registration data to file %s\n',registrationmatfile);
catch exception
  fprintf('Could not save registration data to mat file:\n%s\n',getReport(exception));
end


%
% Save registration_data to text file, for the humans
%
registrationtxtfile = fullfile(expdir,dataloc_params.registrationtxtfilestr);
try
  if exist(registrationtxtfile,'file'),
    delete(registrationtxtfile);
  end
  saveRegistrationDataToTextFile(registrationtxtfile, registration_data) ;  
  fprintf('Saved registration data to txt file %s\n',registrationtxtfile);
catch exception
  fprintf('Could not save registration data to txt file:\n%s\n',getReport(exception));
end

end  % FlyDiscoRegisterTrx()
