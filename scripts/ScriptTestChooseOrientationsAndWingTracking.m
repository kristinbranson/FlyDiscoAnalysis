addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/filehandling;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/simplewing;

intrxfiles = ...
  {'/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/cx_GMR_OL0077B_CsChr_RigD_20150325T155016/fixed_trx1.mat'
  '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/cx_JRC_SS03661_CsChr_RigB_20150325T145528/fixed_trx1.mat'
  '/groups/branson/home/robiea/Projects_data/FlyBubble/OL0046BCsChr_flyBowlAlice_constant_short_20150220T170609/trx.mat'
  '/groups/branson/home/robiea/Projects_data/FlyBubble/OL0046BCsChr_flyBowlAlice_constant_short_20150220T173512/trx.mat'
  '/groups/branson/home/robiea/Projects_data/FlyBubble/OL0046BCsChr_flyBowlAlice_flicker_short_20150220T174440/trx.mat'
  '/groups/branson/home/robiea/Projects_data/FlyBubble/OL0046BCsChr_flyBowlAlice_flicker_short_20150220T175255/trx.mat'
  };
inexpdirs = cellfun(@fileparts,intrxfiles,'Uni',0);

outrootdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';

moviefilestr = 'movie.ufmf';
annfilestr = 'movie.ufmf.ann';
matfilestr = 'trx.mat';
outtrxfilestr = 'testtrx2.mat';
matfilestr2 = 'movie.mat';

paramsfile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/SelectedParamsChooseOrientations20150402.mat';
params = load(paramsfile);

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

%% copy data

outexpdirs = cell(size(inexpdirs));
for i = 1:numel(inexpdirs),
  
  [~,n] = fileparts(inexpdirs{i});
  outexpdirs{i} = fullfile(outrootdir,n);
  
%   if strcmp(outexpdirs{i},inexpdirs{i}),
%     continue;
%   end
  
  if ~exist(outexpdirs{i},'dir'),
    mkdir(outexpdirs{i});
  end
  
  if ~exist(fullfile(outexpdirs{i},moviefilestr),'file'),
    unix(sprintf('ln -s %s/%s %s/%s',inexpdirs{i},moviefilestr,...
      outexpdirs{i},moviefilestr));
  end
  if ~exist(fullfile(outexpdirs{i},annfilestr),'file'),  
    unix(sprintf('ln -s %s/%s %s/%s',inexpdirs{i},annfilestr,...
      outexpdirs{i},annfilestr));
  end
  
  inmatfile = intrxfiles{i};
    
  if exist(fullfile(outexpdirs{i},matfilestr),'file'),
    unix(sprintf('rm %s/%s',outexpdirs{i},matfilestr));
  end
  unix(sprintf('ln -s %s %s/%s',inmatfile,...
    outexpdirs{i},matfilestr));
  
end

%% track these videos with the current settings

firstframe = 1;
endframe = inf;

for expi = 1:numel(outexpdirs),

%   if exist(fullfile(outexpdirs{expi},outtrxfilestr),'file'),
%     continue;
%   end
  
  [trx,pfd,wtinfo,wtunits,trackdata] = ChooseOrientationsAndTrackWings(...
    fullfile(outexpdirs{expi},matfilestr),...
    fullfile(outexpdirs{expi},moviefilestr),...
    'annfile',fullfile(outexpdirs{expi},annfilestr),...
    'paramsfile',paramsfile,...
    'savefile',fullfile(outexpdirs{expi},outtrxfilestr),...
    'debug',0,'firstframe',firstframe,'endframe',endframe,...
    'dofliporientation',false);
  
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
paramsfile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/SelectedParamsChooseOrientations20150701.mat';
save(paramsfile,'-struct','newparams');

%% test

firstframe = 1;
endframe = inf;

for expi = 1:numel(outexpdirs),
  
  tmp = load(fullfile(outexpdirs{expi},outtrxfilestr));
  tmp2 = load(intrxfiles{expi},'trx');
  
  [trx,pfd,wtinfo,wtunits,trackdata] = ChooseOrientationsAndTrackWings(...
    fullfile(outexpdirs{expi},matfilestr),...
    fullfile(outexpdirs{expi},moviefilestr),...
    'annfile',fullfile(outexpdirs{expi},annfilestr),...
    'paramsfile',paramsfile,...
    'savefile','',...
    'debug',0,'firstframe',firstframe,'endframe',endframe,...
    'dofliporientation',true,...
    'trackdata',tmp.trackdata);  

  fprintf('%d, %s: %d\n',expi,outexpdirs{expi},nnz(round(abs([tmp2.trx.theta] - [trx.theta])/pi)==1));

end
