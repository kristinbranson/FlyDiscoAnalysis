FlyDiscodir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis';
analysis_protocol = '20210329_flybubble_FlyBowlRGB_LED';

dataloc_params = ReadParams(fullfile(FlyDiscodir,'settings',analysis_protocol,'dataloc_params.txt'));
% dataloc_params.flytrackertrackstr = 'movie-track.mat';

expdirs_train = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210419_wingtrackingfeatures_testing/socialCsChr_GMR_72C11_AE_01_CsChrimson_RigD_20191114T172654'};

for moviei = 1:numel(expdirs_train),

  expdir = expdirs_train{moviei};
%   [~,expname] = fileparts(expdir);
%   outexpdir = fullfile(rootoutdir,expname);
  FlyTracker2WingTracking(expdir,'dataloc_params',dataloc_params);
  
end