function FlyDiscoMakeCtraxResultsMovie(expdir, varargin)
% Make Ctrax results movie

% Parse the input arguments
[analysis_protocol, settingsdir, datalocparamsfilestr, do_run] = ...
  myparse(varargin,...
          'analysis_protocol','20150915_flybubble_centralcomplex',...
          'settingsdir',default_settings_folder_path(),...
          'datalocparamsfilestr','dataloc_params.txt', ...
          'do_run',[]) ;

% Get locations of parameter files
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

% Get one thing from the registration params
registrationparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.registrationparamsfilestr);
if ~exist(registrationparamsfile,'file'),
  error('Registration params file %s does not exist',registrationparamsfile);
end
raw_registration_params = ReadParams(registrationparamsfile);
registration_params = modernizeRegistrationParams(raw_registration_params) ;
doesYAxisPointUp = registration_params.doesYAxisPointUp ;

% Get one thing from the indicator params
indicatorparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.indicatorparamsfilestr);
if exist(indicatorparamsfile,'file'),
  raw_indicator_params = ReadParams(indicatorparamsfile);
  indicator_params = modernizeIndicatorParams(raw_indicator_params) ;
  isOptogeneticExp = logical(indicator_params.OptogeneticExp) ;
else
  isOptogeneticExp = false ;
end

% Read in the experiment metadata
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
metadata = ReadMetadataFile(metadatafile);

% Determine the ctrax movie parameters file name, and read in the params
defaultctraxresultsmovieparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.ctraxresultsmovieparamsfilestr);
if isOptogeneticExp,
  % Branson Lab optogenetic experiments have a metadata field named
  % "led_protocol", with values like "protocolRGBms_locomotion_GtACR_20220414".
  % The Ctrax results movie params file to use depends on which LED protocol was
  % used.  In particular, the LED protocol files have a date in the file name,
  % and that date matches the date in the "led_protocol" field of the metadata.
  % So we extract that date from the led_protocol field, and check if there's a
  % file named ctraxresultsmovie_params_<date>.txt in the analysis-protocol
  % file.  I think everyone acknowledges this is a little hacky.
  datestrpattern = '20\d{6}';
  match = regexp(metadata.led_protocol,datestrpattern);
  if isempty(match) ,
    % The led_protocol field does not seem to contain a date, so fall back to
    % usual Ctrax results movie parameter file.
    ctraxresultsmovie_params = ReadParams(defaultctraxresultsmovieparamsfile);
    is_using_default_ctrax_results_movie_params = 1 ;
  else
    % The led_protocol field seems to contain a date in the right format, so check
    % for a matching ctraxresultsmovie_params_<date>.txt file.
    ledprotocoldatestr = metadata.led_protocol(match:match+7);
    specificctraxresultsmovie_paramsfile = fullfile(settingsdir,analysis_protocol,['ctraxresultsmovie_params_',ledprotocoldatestr,'.txt']);
    if exist(specificctraxresultsmovie_paramsfile,'file'),
      % The ctraxresultsmovie_params_<date>.txt file exists, so use it.
      ctraxresultsmovie_params = ReadParams(specificctraxresultsmovie_paramsfile);
      is_using_default_ctrax_results_movie_params = 0 ;
    else
      % The ctraxresultsmovie_params_<date>.txt file does not seem to exist, so fall back to
      % usual Ctrax results movie parameter file.
      ctraxresultsmovie_params = ReadParams(defaultctraxresultsmovieparamsfile);
      is_using_default_ctrax_results_movie_params = 1 ;
    end
  end
else
  % If non-optogenetic experiment, use the usual Ctrax results movie parameter file.
  ctraxresultsmovie_params = ReadParams(defaultctraxresultsmovieparamsfile);
  is_using_default_ctrax_results_movie_params = 1 ;
end

% Sort out where the temporary data directory will be
scratch_folder_path = get_scratch_folder_path() ;
defaulttempdatadir = fullfile(scratch_folder_path, 'TempData_FlyDiscoMakeCtraxResultsMovie') ;
if isfield(ctraxresultsmovie_params,'tempdatadir') ,
  tempdatadir = ctraxresultsmovie_params.tempdatadir ;
else
  tempdatadir = defaulttempdatadir ;
end

% Create the temporary data folder if it doesn't exist
if ~exist(tempdatadir, 'dir') ,
  [success1,msg1] = mkdir(tempdatadir);
  if ~success1,
    error('Error making directory %s: %s',tempdatadir,msg1);
  end
end

% Read the trx file if if exists
trxfile = fullfile(expdir,dataloc_params.trxfilestr);
if exist(trxfile, 'file') ,
  trxfilestruct = load(trxfile) ;
  trx = trxfilestruct.trx ;
end

% Read start and end of cropped trajectories
registrationmatfile = fullfile(expdir,dataloc_params.registrationmatfilestr);
if ~exist(registrationmatfile,'file'),
  registration_data = struct() ;
  registration_data.end_frame = max([trx.endframe]);
  registration_data.start_frame = min([trx.firstframe]);
else
  registration_data = load(registrationmatfile);
  if ~isfield(registration_data,'end_frame'),
    registration_data.end_frame = max([trx.endframe]);
  end
  if ~isfield(registration_data,'start_frame'),
    registration_data.start_frame = min([trx.firstframe]);
  end
end

% Read a couple of files if this is an optogenetic experiment
if isOptogeneticExp ,
  indicatorfile = fullfile(expdir,dataloc_params.indicatordatafilestr);
  indicator_data = load(indicatorfile) ;
  ledprotocolfile = fullfile(expdir,dataloc_params.ledprotocolfilestr);
  raw_protocol = loadAnonymous(ledprotocolfile) ;
else
  indicator_data = struct([]) ;
  raw_protocol = struct([]) ;
end

% Make sure the protocol is shorter than the video
error_if_protocol_is_longer_than_video(expdir, settingsdir, analysis_protocol, do_run) ;

% Determine start and end frames of snippets
[firstframes, firstframes_off, endframes_off, nframes, indicatorframes] = ...
  DetermineStartAndEndOfResultsMovieSnippets(isOptogeneticExp, registration_data, ctraxresultsmovie_params, ...
                                             is_using_default_ctrax_results_movie_params, indicator_data, raw_protocol) ;
                                           
% Determine the layout of the results movie frame
[final_nzoomr, final_nzoomc, final_figpos] = ...
  determine_results_movie_figure_layout(ctraxresultsmovie_params, trx) ;
    
% Create subtitle file
[~,basename,~] = fileparts(expdir);
subtitlefile = ...
  write_subtitle_file(expdir, nframes, isOptogeneticExp, ctraxresultsmovie_params, basename, firstframes_off, endframes_off, ...
                      raw_protocol, indicator_data, indicatorframes) ;
                
% Create results movie as .avi
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
avifilestr = sprintf('%s_%s',dataloc_params.ctraxresultsavifilestr,basename);
avi_file_path = fullfile(tempdatadir, [avifilestr, '_temp.avi']);
if exist(avi_file_path,'file'),
  delete(avi_file_path);
end
temp_avi_path = [tempname(scratch_folder_path) '.avi'] ;
[succeeded,~,~,height,width]= ...
  make_ctrax_result_movie('moviename',moviefile,'trxname',trxfile,'aviname',avi_file_path,...
  'nzoomr',final_nzoomr,'nzoomc',final_nzoomc,...
  'boxradius',ctraxresultsmovie_params.boxradius,'taillength',ctraxresultsmovie_params.taillength,...
  'fps',ctraxresultsmovie_params.fps,...
  'maxnframes',nframes,...
  'firstframes',firstframes,...
  'figpos',final_figpos,...
  'movietitle',basename,...
  'compression','none',...
  'useVideoWriter',true,...
  'titletext',false,...
  'showtimestamps',true,...
  'avifileTempDataFile',temp_avi_path,...
  'dynamicflyselection',true,...
  'doshowsex',true, ...
  'doesYAxisPointUp', doesYAxisPointUp);
if ishandle(1),
  close(1);
end
if ~succeeded,
  error('Failed to create Ctrax results movie as .avi file %s',avi_file_path);
end

% Compress the .avi to a .mp4 using the h.264 codec
compressCtraxResultsMovie(expdir, avifilestr, avi_file_path, height, width, subtitlefile) ;
