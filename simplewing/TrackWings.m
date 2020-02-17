function [trx,perframedata,outtrxfile,info] = TrackWings(expdir,varargin)

%% parse parameters

% file names
moviefilestr = 'movie.ufmf';
annfilestr = 'movie.ufmf.ann';
trxfilestr = 'registered_trx.mat';
outtrxfilestr = 'wingtracking_results.mat';
perframedir = 'perframe';
paramsfile = 'WingTrackingParameters.txt';

% first frame to track
firstframe = 1;

% whether to run in debug mode
DEBUG = false;

[moviefilestr,annfilestr,trxfilestr,outtrxfilestr,perframedir,...
  paramsfile,firstframe,DEBUG,isarena,bgmodel,restart,endframe] = ...
  myparse(varargin,...
  'moviefilestr',moviefilestr,...
  'annfilestr',annfilestr,...
  'trxfilestr',trxfilestr,...
  'outtrxfilestr',outtrxfilestr,...
  'perframedir',perframedir,...
  'paramsfile',paramsfile,...
  'firstframe',firstframe,...
  'debug',DEBUG,...
  'isarena',[],...
  'bgmodel',[],...
  'restart',false,...
  'endframe',inf);

if ischar(DEBUG),
  DEBUG = str2double(DEBUG);
end
if ischar(firstframe),
  firstframe = str2double(firstframe);
end

% read parameters from file
params = ReadParams(paramsfile);

%% read stuff
%expdir = fullfile(rootdatadir,experiment_name);

moviefile = fullfile(expdir,moviefilestr);
annfile = fullfile(expdir,annfilestr);
trxfile = fullfile(expdir,trxfilestr);
if isempty(bgmodel),
  [readframe,nframes,fid] = get_readframe_fcn(moviefile); %#ok<ASGLU>
  [nr,nc,~] = size(readframe(1));
  fclose(fid);
  [bgmodel] = read_ann(annfile,'background_center');
  bgmodel = reshape(bgmodel,[nc,nr])';
end
if isempty(isarena),
  [readframe,nframes,fid] = get_readframe_fcn(moviefile); %#ok<ASGLU>
  [nr,nc,~] = size(readframe(1));
  fclose(fid);
  [isarena] = read_ann(annfile,'isarena');
  isarena = reshape(isarena,[nc,nr])' > 0;
end
[trx,~,~,timestamps,trxfiledata] = load_tracks(trxfile); 

tmpfilename = fullfile(expdir,sprintf('TmpResultsTrackWings_%s.mat',datestr(now,'yyyymmddTHHMMSS')));

%% call main function

[trx,perframedata,info,params.units] = TrackWingsHelper(trx,moviefile,bgmodel,isarena,params,...
  'firstframe',firstframe,...
  'debug',DEBUG,...
  'restart',restart,...
  'tmpfilename',tmpfilename,...
  'endframe',endframe);

%% save trx

outtrxfile = fullfile(expdir,outtrxfilestr);
fprintf('Saving trx to file %s...\n',outtrxfile);
timestamp_analyzed = info.trackwings_timestamp; %#ok<NASGU>
version = info.trackwings_version; %#ok<NASGU>
trxfiledata.trx = trx;
trxfiledata.timestamps = timestamps;
trxfiledata.wingtracking_timestamp = timestamp_analyzed;
trxfiledata.wingtracking_version = version;
save(outtrxfile,'-struct','trxfiledata');
%save_tracks(trx,outtrxfile,'timestamps',timestamps);

%% save per-frame data

fprintf('Saving per-frame data...\n');

if ~exist(fullfile(expdir,perframedir),'dir'),
  mkdir(fullfile(expdir,perframedir));
end

fns = fieldnames(perframedata);
for i = 1:numel(fns),
  fn = fns{i};
  s = struct('data',{perframedata.(fn)},'units',params.units.(fn)); %#ok<NASGU>
  filename = fullfile(expdir,perframedir,[fn,'.mat']);
  save(filename,'-struct','s');
end