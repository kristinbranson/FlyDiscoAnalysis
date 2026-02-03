%% setpath
modpath

%% set parameters
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis';
% analysis_protocol = 'settings-internal/20240521_flybubble_LED_VNC3';
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
% analysis_protocol = '20241210_flybubble_LED_VNC3';
% analysis_protocol = '20250721_flybubble_LED_VNC3';
% analysis_protocol = '20251009_flybubble_LED_VNC3'; %same as 20250721_flybubble_LED_VNC3 with updated README
analysis_protocol = '20251009_flybubble_LED_VNC2';
% analysis_protocol = '20251009_flybubble_LED_VNC';

params = {'settingsdir',settingsdir,...
    'analysis_protocol',analysis_protocol, ...
    'forcecompute', false,...
    'debug',true};

%% expdir
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS100086_RigD_20240820T125028'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigB_20240904T124501'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS98080_RigB_20240904T120747'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigB_20240904T124501'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigD_20240529T101336'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS98080_RigB_20240904T120747'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS100086_RigD_20240820T125028'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigB_20240904T124501'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigD_20240529T101336'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_EXT_VGLUT-GAL4_RigC_20240813T120718'};

%test VNC2
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC2_YNA_K_162984_RigB_20230927T120415'};
%failed VNC3
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS97714_RigB_20240919T122459'};
% test VNC
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC_JRC_SS68333_RigB_20210420T151346'};

% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC2_YNA_K_162984_RigB_20230927T120415'};
%failed vglut exps - try to fix
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC2_EXT_VGLUT-GAL4_RigB_20231025T123239'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigB_20240904T124501'};

%reprocess test lines
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS98080_RigB_20240904T120747' ...
%     '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigD_20240529T101336' ...
%     '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_YNA_K_162984_RigB_20240904T124501' ...
%     '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS100086_RigD_20240820T125028'};

% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS100086_RigD_20240820T125028'};
% test vglut empty struct fix
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260202_testinglocostageFIXvglut/VNC2_EXT_VGLUT-GAL4_RigB_20231025T123239'};
explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260112_testing_locomotioncomputeperframestats/VNC2_JRC_SS57983_RigD_20230913T120134'}
%% run code
for i = 1:numel(explist)
    expdir = explist{i};
    FlyDiscoComputeLocomotionMetrics(expdir,params{:})
end
