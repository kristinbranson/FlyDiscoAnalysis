% make results movies
function FlyBowlMakeCtraxResultsMovie(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt');

%% locations of parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% ctrax movie parameters
ctraxresultsmovieparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.ctraxresultsmovieparamsfilestr);
ctraxresultsmovie_params = ReadParams(ctraxresultsmovieparamsfile);
if ~isfield(ctraxresultsmovie_params,'tempdatadir'),
  ctraxresultsmovie_params.tempdatadir = '/groups/branson/bransonlab/projects/olympiad/TempData_FlyBowlMakeCtraxResultsMovie';
end

%% location of data

[~,basename] = fileparts(expdir);
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
trxfile = fullfile(expdir,dataloc_params.trxfilestr);
avifilestr = sprintf('%s_%s',dataloc_params.ctraxresultsavifilestr,basename);
avifile = fullfile(ctraxresultsmovie_params.tempdatadir,[avifilestr,'_temp.avi']);
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

%% create movie

[succeeded,~,~,height,width]= ...
  make_ctrax_result_movie('moviename',moviefile,'trxname',trxfile,'aviname',avifile,...
  'nzoomr',ctraxresultsmovie_params.nzoomr,'nzoomc',ctraxresultsmovie_params.nzoomc,...
  'boxradius',ctraxresultsmovie_params.boxradius,'taillength',ctraxresultsmovie_params.taillength,...
  'fps',ctraxresultsmovie_params.fps,...
  'maxnframes',+ctraxresultsmovie_params.nframes,...
  'firstframes',firstframes,...
  'figpos',ctraxresultsmovie_params.figpos,...
  'movietitle',basename,...
  'compression','none',...
  'useVideoWriter',false,...
  'titletext',false,...
  'avifileTempDataFile',[avifile,'-temp'],...
  'dynamicflyselection',true,...
  'doshowsex',true);

if ishandle(1),
  close(1);
end

if ~succeeded,
  error('Failed to create raw avi %s',avifile);
end

%% create subtitle file

subtitlefile = fullfile(expdir,'subtitles.srt');
fid = fopen(subtitlefile,'w');
dt = [0,ctraxresultsmovie_params.nframes];
ts = cumsum(dt);
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

%% compress

tmpfile = [xvidfile,'.tmp'];
newheight = 4*ceil(height/4);
newwidth = 4*ceil(width/4);
% subtitles are upside down, so encode with subtitles and flip, then flip
% again
cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d,flip -sub %s -subfont-text-scale 2 -msglevel all=2',...
  avifile,tmpfile,newwidth,newheight,subtitlefile);
status = system(cmd);
if status ~= 0,
  fprintf('*****\n');
  warning('Failed to compress avi to %s',xvidfile);
  fprintf('Need to run:\n');
  fprintf('%s\n',cmd);
  cmd2 = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf flip -msglevel all=2',...
    tmpfile,xvidfile);
  fprintf('then\n');
  fprintf('%s\n',cmd2);
  fprintf('then delete %s %s %s\n',tmpfile,avifile,subtitlefile);
  fprintf('*****\n');
else
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
end
