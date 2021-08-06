addpath JAABA/misc;
addpath JAABA/filehandling;
%expdir = '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS46706_RigA_20210512T132556';
rootdir = 'testdata';
% test data during development of code, 20210802
expdir = 'testdata/VNC_JRC_SS46706_RigA_20210512T132556';
analysis_protocol = '20210531_flybubble_LED';
settingsdir = 'settings';

% find an experiment with poor tracking
% vncexpdirs = mydir('/groups/dickson/dicksonlab/flydisco_data','isdir',true);
% nflies = nan(size(vncexpdirs));
% for i = 1:numel(vncexpdirs),
%   vncexpdir = vncexpdirs{i};
%   pff = fullfile(vncexpdir,'perframe','velmag_ctr.mat');
%   if ~exist(pff,'file'),
%     continue;
%   end
%   pd = load(pff);
%   nflies(i) = numel(pd.data);
%   fprintf('%d/%d, nflies = %d\n',i,numel(vncexpdirs),nflies(i));
% end


% test data with poor tracking
% expdir0 = '/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS50831_RigC_20210517T151218';
% %expdir0 = '/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/FTtracked/FlyBubbleRGB/noLED/locomotionGtACR1_29_nonLED_JHS_K_85321_RigA_20210227T113908';
% [~,expname] = fileparts(expdir0);
% expdir = SymbolicCopyExperimentDirectory(expdir0,rootdir);
expdir = 'testdata/VNC_JRC_SS50831_RigC_20210517T151218';

FlyDiscoRegisterTrx(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoDectectIndicatorLedOnOff(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoClassifySex(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoComputePerFrameFeatures(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoComputePerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'docomputehists',true,'debugplot',5);
FlyDiscoPlotPerFrameStats2(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'debug',false);

%%

sd = load(fullfile(expdir,'stats_perframe.mat'));