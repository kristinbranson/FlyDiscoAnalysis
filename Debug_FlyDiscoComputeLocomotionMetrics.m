%% setpath
modpath

%% set parameters
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis';
% analysis_protocol = 'settings-internal/20240521_flybubble_LED_VNC3';
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
% analysis_protocol = '20241210_flybubble_LED_VNC3';
analysis_protocol = '20250721_flybubble_LED_VNC3';
params = {'settingsdir',settingsdir,...
    'analysis_protocol',analysis_protocol, ...
    'forcecompute', true,...
    'debug',true};

%% expdir
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS100086_RigD_20240820T125028'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigB_20240904T124501'};
explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS98080_RigB_20240904T120747'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigB_20240904T124501'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigD_20240529T101336'};
%% run code
for i = 1:numel(explist)
    expdir = explist{i};
    FlyDiscoComputeLocomotionMetrics(expdir,params{:})
end
