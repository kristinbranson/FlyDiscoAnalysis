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

%% location of data

moviefile = fullfile(expdir,dataloc_params.moviefilestr);
trxfile = fullfile(expdir,dataloc_params.trxfilestr);
avifilestr = fullfile(expdir,dataloc_params.ctraxresultsavifilestr);
avifile = [avifilestr,sprintf('_%05dto%05d',[ctraxresultsmovie_params.firstframes;ctraxresultsmovie_params.endframes]),'_temp.avi'];
xvidfile = [avifilestr,sprintf('_%05dto%05d',[ctraxresultsmovie_params.firstframes;ctraxresultsmovie_params.endframes]),'.avi'];

%% create movie

[~,basename] = fileparts(expdir);

succeeded = ...
  make_ctrax_result_movie('moviename',moviefile,'trxname',trxfile,'aviname',avifile,...
  'nzoomr',ctraxresultsmovie_params.nzoomr,'nzoomc',ctraxresultsmovie_params.nzoomc,...
  'boxradius',ctraxresultsmovie_params.boxradius,'taillength',ctraxresultsmovie_params.taillength,...
  'fps',ctraxresultsmovie_params.fps,...
  'maxnframes',ctraxresultsmovie_params.endframes-ctraxresultsmovie_params.firstframes+1,...
  'firstframes',ctraxresultsmovie_params.firstframes,...
  'figpos',ctraxresultsmovie_params.figpos,...
  'movietitle',basename,...
  'compression','none',...
  'useVideoWriter',false,...
  'titletext',false,...
  'avifileTempDataFile',[avifile,'-temp']);

if ishandle(1),
  close(1);
end

if ~succeeded,
  error('Failed to create raw avi %s',avifile);
end

%% create subtitle file

subtitlefile = fullfile(expdir,'subtitles.srt');
fid = fopen(subtitlefile,'w');
dt = [0,ctraxresultsmovie_params.endframes - ctraxresultsmovie_params.firstframes + 1];
ts = cumsum(dt);
for i = 1:numel(t0s),
  fprintf(fid,'%d\n',i);
  fprintf(fid,'%s --> %s\n',...
    datestr(ts(i)/ctraxresultsmovie_params.fps/(3600*24),'HH:MM:SS,FFF'),...
    datestr((ts(i+1)-1)/ctraxresultsmovie_params.fps/(3600*24),'HH:MM:SS,FFF'));
  fprintf(fid,'%s, fr %d-%d\n\n',basename,...
    ctraxresultsmovie_params.firstframes(i),...
    ctraxresultsmovie_params.endframes(i));
end
fclose(fid);

%% compress

cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts bitrate=3000000 -sub %s -subfont-text-scale 2',avifile,xvidfile,subtitlefile);
status = system(cmd);
if status == 0,
  delete(avifile);
  delete(subtitlefile);
else
  error('Failed to compress avi to %s',xvidfile);
end