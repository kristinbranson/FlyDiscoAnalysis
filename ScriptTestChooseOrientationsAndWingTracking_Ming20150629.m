addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/filehandling;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/misc;

outexpdirs = ...
  {'/groups/branson/bransonlab/projects/FlyBowl/data/MingForWingBasedOrientationAlgorithm20150629/GMR_25E04_AE_01_LC16ET3_flyBowlMing_protocolVer3_Ming_LC16epistasis_20150331T125835'
  '/groups/branson/bransonlab/projects/FlyBowl/data/MingForWingBasedOrientationAlgorithm20150629/GMR_82E12_AE_01_LC16ET3_flyBowlMing_protocolVer3_Ming_LC16epistasis_20150331T113930'
  '/groups/branson/bransonlab/projects/FlyBowl/data/MingForWingBasedOrientationAlgorithm20150629/MDN-1LC16ET3_flyBowlMing_protocolVer3_Rai_CsChrimson_2intensities_20150305T153701'
  '/groups/branson/bransonlab/projects/FlyBowl/data/MingForWingBasedOrientationAlgorithm20150629/pBDPGAL4ULC16ET3_flyBowlMing_protocolVer3_Rai_CsChrimson_2intensities_20150305T152132'
  '/groups/branson/bransonlab/projects/FlyBowl/data/MingForWingBasedOrientationAlgorithm20150629/pBDPGAL4U_LC16ET3_flyBowlMing_protocolVer3_Ming_LC16epistasis_20150330T164252'};

moviefilestr = 'movie.ufmf';
annfilestr = 'movie.ufmf.ann';
matfilestr = 'registered_trx.mat';
outtrxfilestr = 'testtrx.mat';

paramsfile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/SelectedParamsChooseOrientations20150402.mat';
params = load(paramsfile);

wingtracking_paramsfile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings/current/WingTrackingParameters.txt';

% params = struct;
% params.min_jump_speed = 5.7256;
% params.nframes_speed = 2;
% params.min_speed = 0.2264;
% params.min_wing_area = 13.0844;
% params.max_wing_area = 22.7860;
% params.max_ecc_confident = 0.6102;
% params.min_ecc_factor = 0.01;
% params.lambda_phi = 0;
% params.lambda_wingarea = 1;
% params.lambda_theta = 10;
% params.wingtracking_params = DefaultWingTrackingParams();
% 
% paramsfile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/ParamsChooseOrientations20150402.mat';
% save(paramsfile,'-struct','params');

%% track these videos with the current settings -- errors have been fixed, so don't allow flipping

firstframe = 1;
endframe = inf;

for expi = 1:numel(outexpdirs),
  
  if exist(fullfile(outexpdirs{expi},outtrxfilestr),'file'),
    continue;
  end
  
  [trx,pfd,wtinfo,wtunits,trackdata] = ChooseOrientationsAndTrackWings(...
    fullfile(outexpdirs{expi},matfilestr),...
    fullfile(outexpdirs{expi},moviefilestr),...
    'paramsfile',paramsfile,...
    'savefile',fullfile(outexpdirs{expi},outtrxfilestr),...
    'debug',1,'firstframe',firstframe,'endframe',endframe,...
    'dofliporientation',false,...
    'wingtracking_params',wingtracking_paramsfile);
  
  save(fullfile(outexpdirs{expi},'trx1.mat'),'trx');
  
  
end

%% combine data

trx = [];
trackdata = {};

trackdatafns = {'wing_areal','wing_arear','istouching'};

for expi = 1:numel(outexpdirs),
  resfile = fullfile(outexpdirs{expi},outtrxfilestr);
  td = load(resfile);
  if expi == 1,
    trx = td.trx;
    trackdata = td.trackdata;
  else
    trx = structappend(trx,td.trx);
    for i = 1:numel(trackdatafns),
      trackdata.(trackdatafns{i}) = [trackdata.(trackdatafns{i}),td.trackdata.(trackdatafns{i})];
    end
  end
end

%% select parameters

newparams = SelectParametersChooseOrientationsAndTrackWings(trx,trackdata);
paramsfile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/SelectedParamsChooseOrientationsMing20150629.mat';
save(paramsfile,'-struct','newparams');

%% test

firstframe = 1;
endframe = inf;

for expi = 1:numel(outexpdirs),
  
  [trx,pfd,wtinfo,wtunits,trackdata] = ChooseOrientationsAndTrackWings(...
    fullfile(outexpdirs{expi},matfilestr),...
    fullfile(outexpdirs{expi},moviefilestr),...
    'paramsfile',paramsfile,...
    'savefile',fullfile(outexpdirs{expi},'testtrx2.mat'),...
    'debug',0,'firstframe',firstframe,'endframe',endframe,...
    'dofliporientation',true,...
    'wingtracking_params',wingtracking_paramsfile);
  
end

for expi = 1:numel(outexpdirs),
  old = load(fullfile(outexpdirs{expi},matfilestr));
  retracked = load(fullfile(outexpdirs{expi},'testtrx2.mat'));
  fprintf('%d, %s: %d\n',expi,outexpdirs{expi},nnz(round(abs([old.trx.theta] - [retracked.trx.theta])/pi)==1));
end