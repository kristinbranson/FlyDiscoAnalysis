% make apt results movie
% uses 3 possible param files:
% aptresultsmovie_params sets headlandmark and taillandmark but can also be input directly. defaults in movie making and 1 and 7 (for alice's fly tracker)
% override_datalocparams - can override the fields of datalocparmas for example trk file name 
% override_ctraxparams - by default the same movie parameters as
% ctraxresultsmovie are used, but these can be overridden with this input
% struct for example fig size of output movie 
% Note additional input options for anonymizing movies not using for pipeline.  

function FlyDiscoMakeAPTResultsMovie(expdir,varargin)

%% Parse the input arguments
[analysis_protocol,...
 settingsdir,...
 datalocparamsfilestr,...
 override_datalocparams,...
 override_ctraxparams,...
 apttrk,...
 hidemovietype,...
 nintervals,...
 hidebackgroundextra,...
 outexpdir,...
 ignoretrxfile,...
 headlandmark,...
 taillandmark,...
 ~,...
 ~,...
 do_run] = ...
  myparse(varargin,...
          'analysis_protocol','20150915_flybubble_centralcomplex',...
          'settingsdir',default_settings_folder_path(),...
          'datalocparamsfilestr','dataloc_params.txt',...
          'override_datalocparams',struct(),...
          'override_ctraxparams',struct(),...
          'apttrk',[],...
          'hidemovietype',false,...
          'nintervals',[],...
          'hidebackgroundextra',.05,...
          'outexpdir','',...
          'ignoretrxfile',false,...
          'headlandmark',nan,...
          'taillandmark',nan,...
          'forcecompute',false,...
          'debug',false,...
          'do_run',[]);

%% locations of parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

fns = fieldnames(override_datalocparams);
for i = 1:numel(fns),
  dataloc_params.(fns{i}) = override_datalocparams.(fns{i});
end

%% read in apt results movie params : headlandmark, taillandmark

datalocparamsfilestr = 'dataloc_params.txt';
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
aptresultsmovieparamsfilestr = dataloc_params.aptresultsmovieparamsfilestr;
% name of parameters file
aptresultsmovieparamsfile = fullfile(settingsdir,analysis_protocol,aptresultsmovieparamsfilestr);
if ~exist(aptresultsmovieparamsfile,'file'),
  error('APT results movie params file %s does not exist',aptresultsmovieparamsfile);
end
% read
aptresultsmovie_params = ReadParams(aptresultsmovieparamsfile);

if isfield(aptresultsmovie_params,'headlandmark')
    headlandmark = aptresultsmovie_params.headlandmark;
end
if isfield(aptresultsmovie_params,'taillandmark')
    taillandmark = aptresultsmovie_params.taillandmark;
end

%% Get one thing from the registration params
% name of parameters file
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
  indicator_params = loadIndicatorParams(indicatorparamsfile) ;
  isOptogeneticExp = logical(indicator_params.OptogeneticExp) ;
else
  isOptogeneticExp = false ;
end

%% location of data

if isempty(outexpdir),
  outexpdir = expdir;
end

[~,basename,~] = fileparts(expdir);
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
if isfield(override_datalocparams,'apttrkfile'),
  apttrkfile = override_datalocparams.apttrkfile;
else
  apttrkfile = fullfile(expdir,dataloc_params.apttrkfilestr);
end
assert(exist(apttrkfile,'file')>0);
if ignoretrxfile,
  trxfile = '';
else
  trxfile = fullfile(expdir,dataloc_params.trxfilestr);
  assert(exist(trxfile,'file')>0);
end
avifilestr = sprintf('%s_%s',dataloc_params.aptresultsavefilestr,basename);
indicatorfile = fullfile(expdir,dataloc_params.indicatordatafilestr);
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
metadata = ReadMetadataFile(metadatafile);
if isOptogeneticExp,
  datestrpattern = '20\d{6}';
  match = regexp(metadata.led_protocol,datestrpattern);
  assert(numel(match)==1);
  
  ledprotocoldatestr = metadata.led_protocol(match:match+7);
  %   ledprotocolfile = fullfile(expdir,'protocol.mat');
  ledprotocolfile = fullfile(expdir,dataloc_params.ledprotocolfilestr);
end

%% ctrax=apt movie parameters
defaultctraxresultsmovieparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.ctraxresultsmovieparamsfilestr);
if isOptogeneticExp,
  specificctraxresultsmovie_paramsfile = fullfile(settingsdir,analysis_protocol,['ctraxresultsmovie_params_',ledprotocoldatestr,'.txt']);
  defaultparams = 1;
  if exist(specificctraxresultsmovie_paramsfile,'file'),
    ctraxresultsmovie_params = ReadParams(specificctraxresultsmovie_paramsfile);
    defaultparams = 0;
  else
    ctraxresultsmovie_params = ReadParams(defaultctraxresultsmovieparamsfile);
  end
else
%  ctraxresultsmovie_params_nonoptogenetic.txt doesn't exist as far as I can tell   
  specificctraxresultsmovie_paramsfile = fullfile(settingsdir,analysis_protocol,'ctraxresultsmovie_params_nonoptogenetic.txt');
  ctraxresultsmovie_params = ReadParams(specificctraxresultsmovie_paramsfile);
  defaultparams = 0;
end


fns = fieldnames(override_ctraxparams);
for i = 1:numel(fns),
  ctraxresultsmovie_params.(fns{i}) = override_ctraxparams.(fns{i});
end

% Sort out where the temporary data directory will be
scratch_folder_path = get_scratch_folder_path() ;
defaulttempdatadir = fullfile(scratch_folder_path, 'TempData_FlyBowlMakeAPTResultsMovie') ;
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

%% read start and end of cropped trajectories

registrationtxtfile = fullfile(expdir,dataloc_params.registrationtxtfilestr);
assert(exist(apttrkfile,'file')>0);
if isempty(apttrk),
  apttrk = TrkFile.load(apttrkfile);
end
if ~exist(registrationtxtfile,'file'),
  registration_params.end_frame = max(apttrk.endframes);
  registration_params.start_frame = min(apttrk.startframes);
else
  registration_params = ReadParams(registrationtxtfile);
  if ~isfield(registration_params,'end_frame'),
    registration_params.end_frame = max(apttrk.endframes);
  end
  if ~isfield(registration_params,'start_frame'),
    registration_params.start_frame = min(apttrk.startframes);
  end
end

%% Make sure the protocol is shorter than the video

error_if_protocol_is_longer_than_video(expdir, settingsdir, analysis_protocol, do_run) ;

%% determine start and end frames of snippets

if ~isOptogeneticExp,
  nframes = registration_params.end_frame-registration_params.start_frame + 1;
  firstframes_off = min(max(0,round(ctraxresultsmovie_params.firstframes*nframes)),nframes-1);
  firstframes_off(ctraxresultsmovie_params.firstframes < 0) = nan;
  middleframes_off = round(ctraxresultsmovie_params.middleframes*nframes);
  middleframes_off(ctraxresultsmovie_params.middleframes < 0) = nan;
  endframes_off = round(ctraxresultsmovie_params.endframes*nframes);
  endframes_off(ctraxresultsmovie_params.endframes < 0) = nan;
  idx = ~isnan(middleframes_off);
  firstframes_off(idx) = ...
   min(nframes-1,max(0,middleframes_off(idx) - ceil(ctraxresultsmovie_params.nframes(idx)/2)));
  idx = ~isnan(endframes_off);
  firstframes_off(idx) = ...
   min(nframes-1,max(0,endframes_off(idx) - ctraxresultsmovie_params.nframes(idx)));
  endframes_off = firstframes_off + ctraxresultsmovie_params.nframes - 1;
  firstframes = registration_params.start_frame + firstframes_off;
else
  if defaultparams,
    load(indicatorfile, 'indicatorLED') ;
    rgbProtocol = loadSingleVariableAnonymously(ledprotocolfile, 'protocol') ;
    protocol = downmixProtocolIfNeeded(rgbProtocol, indicator_params) ;
    
    nsteps = numel(protocol.stepNum);
    % if there's only 1 step 
    if nsteps == 1,
      iter = protocol.iteration;
      n = ceil(iter/6);     % take every nth iteration of the stimulus
      indicatorframes = 1:n:iter;
    elseif nsteps <= 3,      
      % assumes steps have more 2 or more iterations
      if ~all([protocol.iteration] > 1)
          error('Step with iteration < 2. User needs to make a specific ctrax results movie param file')
      end
      indicatorframes = zeros(1,2*nsteps);
      for step = 1:nsteps,
        iter = protocol.iteration(step);
        if step==1,
          indicatorframes(1)=1;
          indicatorframes(2)=iter;
%         this logic doesn't work properly
%         elseif                       
%           indicatorframes(2*step-1)=indicatorframes(2*step-2)+1;
%           indicatorframes(2*step)=indicatorframes(2*step-1);
        elseif step == 2
            indicatorframes(3) = indicatorframes(2)+1;
            indicatorframes(4) = indicatorframes(2)+iter;
        elseif step == 3
            indicatorframes(5) = indicatorframes(4)+1;
            indicatorframes(6) = indicatorframes(4)+iter;

        end        
      end
    else
      n = ceil(nsteps/6);
      index = 1;
      indicatorframes=ones(1,numel(1:n:nsteps));
      for step = 1:n:nsteps,
        if step==1,
          indicatorframes(index)=1;
          index=index+1;
        else
          indicatorframes(index) = indicatorframes(index-1);
          for i = step-n:step-1
            indicatorframes(index) = indicatorframes(index)+protocol.iteration(i);
          end
          index=index+1;
          
        end
      end
    end
    
    % make sure none of the steps are blank (have no stimulus)
    stim = 0;
    step = 0;
    for i=1:numel(indicatorframes),
      while indicatorframes(i) > stim,
        step=step+1;
        stim = stim+protocol.iteration(step);
      end
      if (protocol.intensity(step) == 0),
        if ~(i==numel(indicatorframes)),
          error('Step with intensity = 0 in middle of experiment. User needs to make specific ctrax results movie params file')
        end
        indicatorframes(i) = [];
      end
    end
    
    ctraxresultsmovie_params.indicatorframes = indicatorframes;
    firstframes_off = indicatorLED.startframe(indicatorframes) - ctraxresultsmovie_params.nframes_beforeindicator;  
    % need to do this for specific files as well
    ctraxresultsmovie_params.nframes = ones(1,length(firstframes_off))*ctraxresultsmovie_params.nframes(1);
    endframes_off = firstframes_off + ctraxresultsmovie_params.nframes -1;
    firstframes = registration_params.start_frame + firstframes_off;
  else
    if exist(indicatorfile,'file')
      load(indicatorfile, 'indicatorLED');
      firstframes_off = indicatorLED.startframe(ctraxresultsmovie_params.indicatorframes) - ctraxresultsmovie_params.nframes_beforeindicator;
      endframes_off = firstframes_off + ctraxresultsmovie_params.nframes -1 ;
      firstframes = registration_params.start_frame + firstframes_off;
      ctraxresultsmovie_params.nframes = ones(1,length(firstframes_off))*ctraxresultsmovie_params.nframes;
    end
  end
end
nframesplot = endframes_off - firstframes_off + 1;

if hidemovietype && ~isempty(nintervals) && numel(firstframes) > nintervals,
  intervalidxkeep = round(linspace(1,numel(firstframes),nintervals));
  firstframes = firstframes(intervalidxkeep);
  nframesplot = nframesplot(intervalidxkeep);
end

%% option to not specify nzoomr, nzoomc

if ischar(ctraxresultsmovie_params.nzoomr) || ischar(ctraxresultsmovie_params.nzoomc),
  
  firstframe = min(apttrk.startframes);
  endframe = max(apttrk.endframes);
  trxnframes = endframe-firstframe+1;
  nflies = zeros(1,trxnframes);
  for i = 1:apttrk.ntlts,
    j0 = apttrk.startframes(i)-firstframe+1;
    j1 = apttrk.endframes(i)-firstframe+1;
    nflies(j0:j1) = nflies(j0:j1)+1;
  end
  mediannflies = median(nflies);

  if isnumeric(ctraxresultsmovie_params.nzoomr),
    nzoomr = ctraxresultsmovie_params.nzoomr;
    nzoomc = round(mediannflies/nzoomr);
  elseif isnumeric(ctraxresultsmovie_params.nzoomc),
    nzoomc = ctraxresultsmovie_params.nzoomc;
    nzoomr = round(mediannflies/nzoomc);
  else
    nzoomr = ceil(sqrt(mediannflies));
    nzoomc = round(mediannflies/nzoomr);
  end
  ctraxresultsmovie_params.nzoomr = nzoomr;
  ctraxresultsmovie_params.nzoomc = nzoomc;
  
  if iscell(ctraxresultsmovie_params.figpos),  
    [readframe,~,fid] = get_readframe_fcn(moviefile);
    im = readframe(1);
    [nr,nc,~] = size(im);
    
    rowszoom = floor(nr/nzoomr);
    imsize = [nr,nc+rowszoom*nzoomc];
    figpos = str2double(ctraxresultsmovie_params.figpos);
    if isnan(figpos(3)),
      figpos(3) = figpos(4)*imsize(2)/imsize(1);
    elseif isnan(figpos(4)),
      figpos(4) = figpos(3)*imsize(1)/imsize(2);
    end
    ctraxresultsmovie_params.figpos = figpos;
    
    if fid > 1,
      fclose(fid);
    end
  end
  
end

%% create subtitle file

if ~hidemovietype,
  subtitlefile = fullfile(outexpdir,'subtitles.srt');
  if exist(subtitlefile,'file'),
    delete(subtitlefile);
  end
  fid = fopen(subtitlefile,'w');
  dt = [0,ctraxresultsmovie_params.nframes];
  ts = cumsum(dt);
  
  if ~isOptogeneticExp,
    for i = 1:numel(dt)-1,
      fprintf(fid,'%d\n',i);
      fprintf(fid,'%s --> %s\n',...
        datestr(ts(i)/ctraxresultsmovie_params.fps/(3600*24),'HH:MM:SS,FFF'),...
        datestr((ts(i+1)-1)/ctraxresultsmovie_params.fps/(3600*24),'HH:MM:SS,FFF'));
      fprintf(fid,'%s, fr %d-%d\n\n',basename,...
        firstframes_off(i)+1,...
        endframes_off(i)+1);
    end
    fclose(fid);
  else
    
    protocolOfSomeKind = loadSingleVariableAnonymously(ledprotocolfile, 'protocol');
    protocol = downmixProtocolIfNeeded(protocolOfSomeKind, indicator_params) ;
    
    load(indicatorfile, 'indicatorLED') ;
    stimtimes = indicatorLED.starttimes(ctraxresultsmovie_params.indicatorframes);
    
    j = 1;
    t = protocol.duration(j)/1000;
    
    step = zeros(1,length(stimtimes));
    freq = zeros(1,length(stimtimes));
    intensity = zeros(1,length(stimtimes));
    dutycycle = zeros(1,length(stimtimes));
    duration = zeros(1,length(stimtimes));
    
    for i = 1:length(stimtimes),
      while stimtimes(i) > t
        j = j+1;
        t = t + protocol.duration(j)/1000;
      end
      
      step(i) = j;
      % not sure this logic is good
      if protocol.pulseNum(j) > 1,
        freq(i) = 1/(protocol.pulsePeriodSP(j)/1000);
        stim_type{i} = ['Plsd ', num2str(freq(i)), 'Hz']; %#ok<AGROW> 
      else
        stim_type{i} = 'Cnst'; %#ok<AGROW> 
      end
      
      intensity(i) = protocol.intensity(j);
      dutycycle(i) = (protocol.pulseWidthSP(j)/protocol.pulsePeriodSP(j))*100;
      duration(i) = protocol.pulseNum(j)*protocol.pulsePeriodSP(j);
      
    end
    
    for k = 1:numel(dt)-1,
      fprintf(fid,'%d\n',k);
      
      t_start = ts(k)/ctraxresultsmovie_params.fps/(3600*24);
      t_end = (ts(k) + min(ctraxresultsmovie_params.fps,(ts(k+1)-1-ts(k))/ctraxresultsmovie_params.subdecimationfactor))/ ...
        ctraxresultsmovie_params.fps/(3600*24);
      
      fprintf(fid,'%s --> %s\n',...
        datestr(t_start,'HH:MM:SS,FFF'),...
        datestr(t_end,'HH:MM:SS,FFF'));
      fprintf(fid,'%s\n%s %s %s %s %s %s\n\n',basename,...
        ['Step ', num2str(step(k))],...
        ['(Stim ', num2str(ctraxresultsmovie_params.indicatorframes(k)),'/', ...
        num2str(numel(indicatorLED.starttimes)),'):'],...
        stim_type{k},...
        ['(',num2str(dutycycle(k)),'% on)'],...
        ['at ',num2str(intensity(k)), '% intensity'],...
        ['for ',num2str(duration(k)), ' ms']);
      
    end
    fclose(fid);
  end
else
  subtitlefile = '';
end

%% mask to hide identifying information from movie

if hidemovietype,
  
  rd = load(fullfile(expdir,dataloc_params.registrationmatfilestr));
  [readframe,~,fid] = get_readframe_fcn(moviefile);
  im = readframe(1);
  [nr,nc,~] = size(im);
  if fid > 1,
    fclose(fid);
  end

  [xgrid,ygrid] = meshgrid(1:nc,1:nr);
  dcenter = sqrt((xgrid-rd.circleCenterX).^2 + (ygrid-rd.circleCenterY).^2)./rd.circleRadius;
  dohide = dcenter >= 1+hidebackgroundextra;
else
  dohide = [];  
end



%% create movie
avi_file_path = fullfile(tempdatadir, [avifilestr, '_temp.avi']);
if exist(avi_file_path,'file'),
  delete(avi_file_path);
end
temp_avi_path = [tempname(scratch_folder_path) '.avi'] ;
[succeeded,~,~,height,width]= ...
  make_apt_result_movie('moviename',moviefile,'aptname',apttrk,'trxname',trxfile,'aviname',avi_file_path,...
  'nzoomr',ctraxresultsmovie_params.nzoomr,'nzoomc',ctraxresultsmovie_params.nzoomc,...
  'boxradius',ctraxresultsmovie_params.boxradius,'taillength',ctraxresultsmovie_params.taillength,...
  'fps',ctraxresultsmovie_params.fps,...
  'maxnframes',nframesplot,...
  'firstframes',firstframes,...
  'figpos',ctraxresultsmovie_params.figpos,...
  'movietitle',basename,...
  'compression','none',...
  'useVideoWriter',true,...
  'titletext',false,...
  'showtimestamps',true,...
  'avifileTempDataFile',temp_avi_path,...
  'dynamicflyselection',true,...
  'doshowsex',~hidemovietype,...
  'hidemovietype',hidemovietype,...
  'hidemask',dohide,...
  'hidemaskvalue',0,...
  'headlandmark',headlandmark,...
  'taillandmark',taillandmark, ...
  'doesYAxisPointUp', doesYAxisPointUp);
%'fps',ctraxresultsmovie_params.fps,...
    %'maxnframes',+ctraxresultsmovie_params.nframes,...
if ishandle(1),
  close(1);
end

if ~succeeded,
  error('Failed to create raw avi %s',avi_file_path);
end

%% compress

%tmpfile = [xvidfile,'.tmp'];
newheight = 4*ceil(height/4);
newwidth = 4*ceil(width/4);


% subtitles are upside down, so encode with subtitles and flip, then flip
% again
%cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d,flip -sub %s -subfont-text-scale 2 -msglevel all=2',...
%  avifile,tmpfile,newwidth,newheight,subtitlefile);

mp4_file_path = fullfile(outexpdir, [avifilestr,'.mp4']) ;

nowstr = datestr(now,'yyyymmddTHHMMSSFFF');
passlogfile = sprintf('%s_%s',avi_file_path,nowstr);
ffmpeg_command = 'env -u LD_LIBRARY_PATH /usr/bin/ffmpeg' ;  
    % Use env -u to clear Matlab's very Matlab-specific
    % LD_LIBRARY_PATH
% if isequal(get_distro_codename(), 'Ubuntu') && exist('/usr/bin/ffmpeg', 'file') ,
%     ffmpeg_command = 'env -u LD_LIBRARY_PATH /usr/bin/ffmpeg' ;
%       % Use the local ffmpeg, to avoid fontconfig issues
%       % Have to use env -u to clear Matlab's very Matlab-specific
%       % LD_LIBRARY_PATH
% else
%     %ffmpeg_command = '/misc/local/ffmpeg-4.3.1/bin/ffmpeg' ;  % cluster pre OL9
%     ffmpeg_command = '/misc/sc/ffmpeg-git-20230313-amd64-static/ffmpeg' ;
% end
if hidemovietype,
  subtitlestr = '';
else
  subtitlestr = sprintf(' -vf "subtitles=%s:force_style=''FontSize=10,FontName=Helvetica''"',subtitlefile);
end
cmd = sprintf('%s -i %s -y -passlogfile %s -c:v h264 -pix_fmt yuv420p -s %dx%d -b:v 1600k%s -pass 1 -f mp4 /dev/null',...
  ffmpeg_command, avi_file_path,passlogfile,newwidth,newheight,subtitlestr);
cmd2 = sprintf('%s -i %s -y -passlogfile %s -c:v h264 -pix_fmt yuv420p -s %dx%d -b:v 1600k%s -pass 2 -f mp4 %s',...
  ffmpeg_command, avi_file_path,passlogfile,newwidth,newheight,subtitlestr,mp4_file_path);

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
    if ~hidemovietype, delete(subtitlefile); end
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