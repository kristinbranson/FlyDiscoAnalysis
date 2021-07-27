addpath JAABA/misc;
addpath JAABA/filehandling;
%expdir = '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS46706_RigA_20210512T132556';
expdir = 'testdata/VNC_JRC_SS46706_RigA_20210512T132556';
analysis_protocol = '20210531_flybubble_LED';
settingsdir = 'settings';

% FlyDiscoRegisterTrx(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoDectectIndicatorLedOnOff(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoClassifySex(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoComputePerFrameFeatures(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoComputePerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'docomputehists',true,'debugplot',0);

%%

sd = load(fullfile(expdir,'stats_perframe.mat'));