FlyDiscodir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis';
analysis_protocol = '20210329_flybubble_FlyBowlRGB_LED';

dataloc_params = ReadParams(fullfile(FlyDiscodir,'settings',analysis_protocol,'dataloc_params.txt'));
% dataloc_params.flytrackertrackstr = 'movie-track.mat';

expdirs_train = {'/groups/branson/home/robiea/Projects_data/FlyDisco/KatieTestData/FlyBowlDisco_RGBonly_401/20210401T132850_rig1_flyBowl1__aIPgSS1UASCsChrimson_KS_redonly_protocolRGB_0315_2'};

for moviei = 1:numel(expdirs_train),

  expdir = expdirs_train{moviei};
%   [~,expname] = fileparts(expdir);
%   outexpdir = fullfile(rootoutdir,expname);
  FlyTracker2WingTracking(expdir,'dataloc_params',dataloc_params);
  
end