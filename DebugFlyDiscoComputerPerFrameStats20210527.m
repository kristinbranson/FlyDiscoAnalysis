addpath JAABA/misc;
addpath JAABA/filehandling;
%expdir = '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS46706_RigA_20210512T132556';
rootdir = 'testdata';
% test data during development of code, 20210802
expdir = 'testdata/VNC_JRC_SS46706_RigA_20210512T132556';
analysis_protocol = '20210531_flybubble_LED';
analysis_protocol = '20210531_flybubble_LED_AR_20210819';
settingsdir = 'settings';
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings';

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

% test data selected by Ryo
%expdir0 = '/groups/branson/bransonlab/flydisco_data/VNC_EXT_VGLUT-GAL4_RigA_20210427T125905';
%expdir0 = '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS68333_RigA_20210422T150926';
expdir0 = '/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS62014_RigD_20210525T133656';
expdir0 = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210820_testingperframe/VNC_EXT_VGLUT-GAL4_RigC_20210628T152511';

[~,expname] = fileparts(expdir0);
expdir = SymbolicCopyExperimentDirectory(expdir0,rootdir);
% results copied to 
% /misc/public/Kristin2Disco/VNC_EXT_VGLUT-GAL4_RigA_20210427T125905/stats.html
% /misc/public/Kristin2Disco/VNC_JRC_SS68333_RigA_20210422T150926/stats.html
% /misc/public/Kristin2Disco/VNC_JRC_SS62014_RigD_20210525T133656/stats.html


FlyDiscoRegisterTrx(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoDetectIndicatorLedOnOff(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoClassifySex(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoComputePerFrameFeatures(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoComputePerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'docomputehists',true,'debugplot',5,'dorecompute',true);
FlyDiscoPlotPerFrameStats2(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'debug',false,...
  'plotstats',2,'plotstim',2,'makestimvideos',2,'plothist',2,'plotflies',true,'plotstimtrajs',2);

%%

sd = load(fullfile(expdir,'stats_perframe.mat'));