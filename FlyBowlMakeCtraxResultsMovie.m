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

%% 

[~,basename] = fileparts(expdir);

succeeded = ...
  make_ctrax_result_movie('moviename',moviefile,'trxname',trxfile,'aviname',avifile,...
  'nzoomr',ctraxresultsmovie_params.nzoomr,'nzoomc',ctraxresultsmovie_params.nzoomc,...
  'boxradius',ctraxresultsmovie_params.boxradius,'taillength',ctraxresultsmovie_params.taillength,...
  'fps',ctraxresultsmovie_params.fps,...
  'maxnframes',ctraxresultsmovie_params.endframes-ctraxresultsmovie_params.firstframes+1,...
  'firstframes',ctraxresultsmovie_params.firstframes,...
  'figpos',ctraxresultsmovie_params.figpos,...
  'movietitle',basename);

if ~succeeded,
  error('Failed to create raw avi %s',avifile);
end

%%

cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts bitrate=3000000',avifile,xvidfile);
status = system(cmd);
if status == 0,
  delete(avifile);
else
  error('Failed to compress avi to %s',xvidfile);
end