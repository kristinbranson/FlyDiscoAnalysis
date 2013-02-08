% make results movies
function FlyBowlMakeJAABAResultsMovie(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt');
subtitles = false;

%% locations of parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% jaaba movie parameters
jaabaresultsmovieparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.jaabaresultsmovieparamsfilestr);
jaabaresultsmovie_params = ReadParams(jaabaresultsmovieparamsfile);
if ~isfield(jaabaresultsmovie_params,'tempdatadir'),
  jaabaresultsmovie_params.tempdatadir = '/groups/branson/bransonlab/projects/olympiad/TempData_FlyBowlMakeCtraxResultsMovie';
end
classifierparamsfiles = cellfun(@(x) fullfile(settingsdir,analysis_protocol,x),dataloc_params.jaabaclassifierparamsfilestrs,...
  'UniformOutput',false);

if isfield(jaabaresultsmovie_params,'behavior2color'),
  n = numel(jaabaresultsmovie_params.behavior2color);
  behavior2color = cell(floor(n/2),2);
  for i = 1:floor(n/2),
    behavior2color{i,1} = jaabaresultsmovie_params.behavior2color{2*i-1};
    behavior2color{i,2} = str2double(regexp(jaabaresultsmovie_params.behavior2color{2*i},'\s+','split'));
  end
else
  behavior2color = {};
end

%% location of data

[~,basename] = fileparts(expdir);
%moviefile = fullfile(expdir,dataloc_params.moviefilestr);
trxfile = fullfile(expdir,dataloc_params.trxfilestr);
avifilestr = sprintf('%s_%s',dataloc_params.jaabaresultsavifilestr,basename);
avifile = fullfile(jaabaresultsmovie_params.tempdatadir,[avifilestr,'_temp.avi']);
xvidfile = fullfile(expdir,[avifilestr,'.avi']);

%% read start and end of cropped trajectories

registrationtxtfile = fullfile(expdir,dataloc_params.registrationtxtfilestr);
registration_params = ReadParams(registrationtxtfile);
if ~isfield(registration_params,'end_frame'),
  load(trxfile,'trx');
  registration_params.end_frame = max([trx.endframe]);
end
if ~isfield(registration_params,'start_frame'),
  if ~exist('trx','var'),
    load(trxfile,'trx');
  end
  registration_params.start_frame = min([trx.firstframe]);
end
nframes = registration_params.end_frame-registration_params.start_frame + 1;
firstframes_off = min(max(0,round(jaabaresultsmovie_params.firstframes*nframes)),nframes-1);
firstframes_off(jaabaresultsmovie_params.firstframes < 0) = nan;
middleframes_off = round(jaabaresultsmovie_params.middleframes*nframes);
middleframes_off(jaabaresultsmovie_params.middleframes < 0) = nan;
endframes_off = round(jaabaresultsmovie_params.endframes*nframes);
endframes_off(jaabaresultsmovie_params.endframes < 0) = nan;
idx = ~isnan(middleframes_off);
firstframes_off(idx) = ...
  min(nframes-1,max(0,middleframes_off(idx) - ceil(jaabaresultsmovie_params.nframes(idx)/2)));
idx = ~isnan(endframes_off);
firstframes_off(idx) = ...
  min(nframes-1,max(0,endframes_off(idx) - jaabaresultsmovie_params.nframes(idx)));
endframes_off = firstframes_off + jaabaresultsmovie_params.nframes - 1;


firstframes = registration_params.start_frame + firstframes_off;

%% create movie

[succeeded,~,~,height,width]= ...
  make_jaaba_result_movie_plotallflies(expdir,'moviefilestr',dataloc_params.moviefilestr,'trxfilestr',dataloc_params.trxfilestr,'aviname',avifile,...
  'nzoomr',jaabaresultsmovie_params.nzoomr,'nzoomc',jaabaresultsmovie_params.nzoomc,...
  'boxradius',jaabaresultsmovie_params.boxradius,'taillength',jaabaresultsmovie_params.taillength,...
  'fps',jaabaresultsmovie_params.fps,...
  'maxnframes',+jaabaresultsmovie_params.nframes,...
  'firstframes',firstframes,...
  'figpos',jaabaresultsmovie_params.figpos,...
  'movietitle',basename,...
  'compression','none',...
  'useVideoWriter',false,...
  'titletext',true,...
  'avifileTempDataFile',[avifile,'-temp'],...
  'dynamicflyselection',true,...
  'doshowsex',true,...
  'classifierparamsfiles',classifierparamsfiles,...
  'behavior2color',behavior2color,...
  'debug',false);

if ishandle(1),
  close(1);
end

if ~succeeded,
  error('Failed to create raw avi %s',avifile);
end

%% create subtitle file

if subtitles,

subtitlefile = fullfile(expdir,'subtitles.srt');
fid = fopen(subtitlefile,'w');
dt = [0,jaabaresultsmovie_params.nframes];
ts = cumsum(dt);
for i = 1:numel(dt)-1,
  fprintf(fid,'%d\n',i);
  fprintf(fid,'%s --> %s\n',...
    datestr(ts(i)/jaabaresultsmovie_params.fps/(3600*24),'HH:MM:SS,FFF'),...
    datestr((ts(i+1)-1)/jaabaresultsmovie_params.fps/(3600*24),'HH:MM:SS,FFF'));
  fprintf(fid,'%s, fr %d-%d\n\n',basename,...
    firstframes_off(i)+1,...
    endframes_off(i)+1);
end
fclose(fid);

end

%% compress

tmpfile = [xvidfile,'.tmp'];
newheight = 4*ceil(height/4);
newwidth = 4*ceil(width/4);
% subtitles are upside down, so encode with subtitles and flip, then flip
% again
if subtitles,
  cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d,flip -sub %s -subfont-text-scale 2 -msglevel all=2',...
    avifile,tmpfile,newwidth,newheight,subtitlefile);
else
  cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d -msglevel all=2',...
    avifile,xvidfile,newwidth,newheight);
end

status = system(cmd);
if status ~= 0,
  fprintf('*****\n');
  warning('Failed to compress avi to %s',xvidfile);
  fprintf('Need to run:\n');
  fprintf('%s\n',cmd);
  if subtitles,
    cmd2 = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf flip -msglevel all=2',...
      tmpfile,xvidfile);
    fprintf('then\n');
    fprintf('%s\n',cmd2);
    fprintf('then delete %s %s %s\n',tmpfile,avifile,subtitlefile);
  else
    fprintf('then delete %s\n',avifile);
  end

  fprintf('*****\n');
else
  if subtitles,
    cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf flip -msglevel all=2',...
      tmpfile,xvidfile);
    status = system(cmd);
    if status ~= 0,
      fprintf('*****\n');
      warning('Failed to add subtitles to %s',xvidfile);
      fprintf('Need to run:\n');
      fprintf('%s\n',cmd);
      fprintf('then delete %s %s %s\n',tmpfile,avifile,subtitlefile);
      fprintf('*****\n');
    else
      delete(tmpfile);
      delete(avifile);
      delete(subtitlefile);
    end
  else
    delete(avifile);
  end
end
