function FlyDiscoMakeCtraxResultsMovie(expdir,varargin)
% make results movies

[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','20150915_flybubble_centralcomplex',...
  'settingsdir',internal_settings_folder_path(),...
  'datalocparamsfilestr','dataloc_params.txt');

%% locations of parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% location of data

[~,basename,~] = fileparts(expdir);
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
trxfile = fullfile(expdir,dataloc_params.trxfilestr);
avifilestr = sprintf('%s_%s',dataloc_params.ctraxresultsavifilestr,basename);
%xvidfile = fullfile(expdir,[avifilestr,'.avi']);
% metadatafile = fullfile(expdir,'Metadata.xml');
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
metadata = ReadMetadataFile(metadatafile);
% commonregistrationparamsfile = fullfile(settingsdir,analysis_protocol,'registration_params.txt');
commonregistrationparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.registrationparamsfilestr);
commonregistrationparams = ReadParams(commonregistrationparamsfile);
if commonregistrationparams.OptogeneticExp,
  indicatorfile = fullfile(expdir,dataloc_params.indicatordatafilestr);
  datestrpattern = '20\d{6}';
  match = regexp(metadata.led_protocol,datestrpattern);  
  ledprotocoldatestr = metadata.led_protocol(match:match+7);
  ledprotocolfile = fullfile(expdir,dataloc_params.ledprotocolfilestr);
else
  % non-optogenetic experiments don't have these things
  indicatorfile = [] ;
  ledprotocolfile = [] ;  
end

%% ctrax movie parameters
defaultctraxresultsmovieparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.ctraxresultsmovieparamsfilestr);
if commonregistrationparams.OptogeneticExp,
  specificctraxresultsmovie_paramsfile = fullfile(settingsdir,analysis_protocol,['ctraxresultsmovie_params_',ledprotocoldatestr,'.txt']);
  if exist(specificctraxresultsmovie_paramsfile,'file'),
    ctraxresultsmovie_params = ReadParams(specificctraxresultsmovie_paramsfile);
    is_using_default_ctrax_results_movie_params = 0 ;
  else
    ctraxresultsmovie_params = ReadParams(defaultctraxresultsmovieparamsfile);
    is_using_default_ctrax_results_movie_params = 1 ;
  end
else
  %  ctraxresultsmovie_params_nonoptogenetic.txt doesn't exist as far as I can tell
  specificctraxresultsmovie_paramsfile = fullfile(settingsdir,analysis_protocol,'ctraxresultsmovie_params_nonoptogenetic.txt');
  ctraxresultsmovie_params = ReadParams(specificctraxresultsmovie_paramsfile);
  is_using_default_ctrax_results_movie_params = nan ;  % want an error if we ever try to use this in an if statement
end

% Sort out where the temporary data directory will be
scratch_folder_path = get_scratch_folder_path() ;
defaulttempdatadir = fullfile(scratch_folder_path, 'TempData_FlyBowlMakeCtraxResultsMovie') ;
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
if exist(trxfile, 'file') ,
  trxfilestruct = load(trxfile) ;
  trx = trxfilestruct.trx ;
end

%% read start and end of cropped trajectories

registrationtxtfile = fullfile(expdir,dataloc_params.registrationtxtfilestr);
if ~exist(registrationtxtfile,'file'),
  %load(trxfile,'trx');
  registration_params.end_frame = max([trx.endframe]);
  registration_params.start_frame = min([trx.firstframe]);
else
  registration_params = ReadParams(registrationtxtfile);
  if ~isfield(registration_params,'end_frame'),
    %load(trxfile,'trx');
    registration_params.end_frame = max([trx.endframe]);
  end
  if ~isfield(registration_params,'start_frame'),
%     if ~exist('trx','var'),
%       load(trxfile,'trx');
%     end
    registration_params.start_frame = min([trx.firstframe]);
  end
end

% Read a couple of files if this is an optogenetic experiment
if commonregistrationparams.OptogeneticExp ,
  indicatorstruct = load(indicatorfile) ;
else
  indicatorstruct = struct([]) ;
end
if commonregistrationparams.OptogeneticExp ,
  ledprotocolstruct = load(ledprotocolfile) ;
else
  ledprotocolstruct = struct([]) ;
end



% determine start and end frames of snippets
[firstframes, firstframes_off, endframes_off, nframes, indicatorframes] = ...
    DetermineStartAndEndOfResultsMovieSnippets(commonregistrationparams, registration_params, ctraxresultsmovie_params, ...
                                               is_using_default_ctrax_results_movie_params, indicatorstruct, ledprotocolstruct, metadata) ;

                                           
% Determine the layout of the results movie frame
[final_nzoomr, final_nzoomc, final_figpos] = ...
    determine_results_movie_figure_layout(ctraxresultsmovie_params, trx) ;

    
    
% create subtitle file
subtitlefile = ...
    write_subtitle_file(expdir, nframes, commonregistrationparams, ctraxresultsmovie_params, basename, firstframes_off, endframes_off, ...
                        metadata, ledprotocolstruct, indicatorstruct, indicatorframes) ;


                
% Create results movie
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
  'doshowsex',true);
%'fps',ctraxresultsmovie_params.fps,...
    %'maxnframes',+ctraxresultsmovie_params_nframes,...
if ishandle(1),
  close(1);
end

if ~succeeded,
  error('Failed to create raw avi %s',avi_file_path);
end

% compress

%tmpfile = [xvidfile,'.tmp'];
newheight = 4*ceil(height/4);
newwidth = 4*ceil(width/4);


% subtitles are upside down, so encode with subtitles and flip, then flip
% again
%cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d,flip -sub %s -subfont-text-scale 2 -msglevel all=2',...
%  avifile,tmpfile,newwidth,newheight,subtitlefile);

mp4_file_path = fullfile(expdir, [avifilestr,'.mp4']) ;
nowstr = datestr(now,'yyyymmddTHHMMSSFFF');
passlogfile = sprintf('%s_%s',avi_file_path,nowstr);
if isequal(get_distro_codename(), 'Ubuntu') && exist('/usr/bin/ffmpeg', 'file') ,
    ffmpeg_command = 'env -u LD_LIBRARY_PATH /usr/bin/ffmpeg' ;  
        % Use the local ffmpeg, to avoid fontconfig issues
        % Have to use env -u to clear Matlab's very Matlab-specific
        % LD_LIBRARY_PATH
else
    ffmpeg_command = '/misc/local/ffmpeg-4.3.1/bin/ffmpeg' ;
end
cmd = sprintf('%s -i %s -y -passlogfile %s -c:v h264 -pix_fmt yuv420p -s %dx%d -b:v 1600k -vf "subtitles=%s:force_style=''FontSize=10,FontName=Helvetica''" -pass 1 -f mp4 /dev/null',...
  ffmpeg_command, avi_file_path,passlogfile,newwidth,newheight,subtitlefile);
cmd2 = sprintf('%s -i %s -y -passlogfile %s -c:v h264 -pix_fmt yuv420p -s %dx%d -b:v 1600k -vf "subtitles=%s:force_style=''FontSize=10,FontName=Helvetica''" -pass 2 -f mp4 %s',...
  ffmpeg_command, avi_file_path,passlogfile,newwidth,newheight,subtitlefile,mp4_file_path);

status = system(cmd);
if status ~= 0,
  fprintf('*****\n');
  warning('ffmpeg first pass failed.');
  fprintf('Need to run:\n');
  fprintf('%s\n',cmd);
%  cmd2 = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf flip -msglevel all=2',...
%    tmpfile,xvidfile);
  fprintf('then\n');
  fprintf('%s\n',cmd2);
  fprintf('then delete %s %s* %s\n',avi_file_path,passlogfile,subtitlefile);
  fprintf('*****\n');
else
%   cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf flip -msglevel all=2',...
%     tmpfile,xvidfile);
  status = system(cmd2);
  if status ~= 0,
    fprintf('*****\n');
    warning('ffmpeg second pass failed.');
    fprintf('Need to run:\n');
    fprintf('%s\n',cmd2);
    fprintf('then delete %s %s* %s\n',avi_file_path,passlogfile,subtitlefile);
    fprintf('*****\n');    
  else
    %delete(tmpfile);
    delete(avi_file_path);
    delete(subtitlefile);
    fname = [passlogfile '-*.log'];
    if isscalar(dir(fname))
      delete(fname);
    end    
    fname = [passlogfile '-*.log.mbtree'];
    if isscalar(dir(fname))
      delete(fname);
    end    
  end
end