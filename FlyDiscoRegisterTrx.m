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
registration_params = ReadParams(registrationparamsfile);

% Determine final dotemporalreg, dotemporaltruncation
[dotemporalreg, dotemporaltruncation, idealVideoDuration] = determineTemporalStuff(registration_params, dotemporalreg, dotemporaltruncation) ;
if dotemporalreg && dotemporaltruncation ,
  error('Temporal registration and temporal truncation are mutually exclusive.  You can''t do both!') ;
end
if dotemporaltruncation && isempty(idealVideoDuration) ,
  error('Temporal truncation requested, but the idealVideoDuration is []') ;
end


%
% Determine the registration transform.  Along the way, detect the
% registration mark(s), even though they are currently not used for
% registration.
%
registration_data_0 = determineRegistrationTransform(expdir, analysis_protocol_folder_path, dataloc_params, registration_params) ;
fprintf('Determined registration transform.\n');


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
if ~isfield(registration_params,'maxFlyTrackerNanInterpFrames'),
  fprintf('maxFlyTrackerNanInterpFrames not set in registration_params, using default value.\n');
  args = {};
else
  args = {'maxFlyTrackerNanInterpFrames',registration_params.maxFlyTrackerNanInterpFrames};
end
[trx,ninterpframes,newid2oldid] = PostprocessFlyTracker(trx,args{:});
fprintf('Removed nans from tracker output.\n');
fprintf('Number of nans interpolated through: %d frames\n',ninterpframes);
fprintf('Number of identities was %d, now %d\n',nids0,numel(trx));

% Store some metadata about nan-removal in registration_data
registration_data_1 = setfields(registration_data_0, 'flytracker_nnanframes', ninterpframes, 'flytracker_nids0', nids0) ;

% Determine if timestamps are reliable, and compute a fallback dt if not
[are_timestamps_reliable, fallback_dt] = getFallbackDtIfNeeded(registration_params, timestamps, trx) ;

% If there are zero tracks, something is wrong
if isempty(trx),
  error('No flies tracked.');
end

% Apply spatial registration to trajectories
trx = appendPhysicalUnitFieldsToTrx(trx, registration_data_1, are_timestamps_reliable, fallback_dt) ;
fprintf('Applied spatial registration.\n') ;



%
% Crop start and end of trajectories based on fly loaded time
%
[trx, timestamps, registration_data_2, newid2oldid] = ...
  performTemporalRegistration(dotemporalreg, trx, timestamps, registration_data_1, newid2oldid, ...
                              trx_file_names_that_are_not_per_frame, expdir, dataloc_params) ;


% 
% Truncate end of movie based on value from registration params 
%
[trx, registration_data] = ...
  performTemporalTruncation(dotemporaltruncation, trx, registration_data_2, newid2oldid, ...
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
