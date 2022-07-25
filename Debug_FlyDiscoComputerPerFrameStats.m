modpath


% rootdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210820_testingperframe';
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210820_testingperframe/VNC_EXT_VGLUT-GAL4_RigC_20210628T152511'};
% VNC_JRC_SS50051_RigD_20210426T151519
% explist = {'/groups/branson/bransonlab/alice/20211124_testingperframestats/VNC_EXT_VGLUT-GAL4_RigA_20210427T125905',...
% '/groups/branson/bransonlab/alice/20211124_testingperframestats/VNC_JRC_SS44225_RigA_20210928T134724',...
% '/groups/branson/bransonlab/alice/20211124_testingperframestats/VNC_JRC_SS46233_RigA_20210415T145311',...
% '/groups/branson/bransonlab/alice/20211124_testingperframestats/VNC_JRC_SS62014_RigD_20210525T133656',...
% '/groups/branson/bransonlab/alice/20211124_testingperframestats/VNC_JRC_SS68333_RigA_20210422T150926',...
% '/groups/branson/bransonlab/alice/20211124_testingperframestats/VNC_JRC_SS71988_RigA_20210914T143410',...
% '/groups/branson/bransonlab/alice/20211124_testingperframestats/VNC_YNA_K_162984_RigC_20210526T155035'};
% explist = {'/groups/branson/bransonlab/alice/20211124_testingperframestats/VNC_JRC_SS71988_RigA_20210914T143410'};
explist = {'/groups/branson/bransonlab/alice/20211124_testingperframestats/VNC_EXT_VGLUT-GAL4_RigA_20210427T125905'};

% analysis_protocol = '20210531_flybubble_LED';
% analysis_protocol = '20210531_flybubble_LED_AR_20210819';
analysis_protocol = '20210531_flybubble_LED_AR_20220405';
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
% expdir = 'testdata/VNC_JRC_SS50831_RigC_20210517T151218';

% test data selected by Ryo
%expdir0 = '/groups/branson/bransonlab/flydisco_data/VNC_EXT_VGLUT-GAL4_RigA_20210427T125905';
%expdir0 = '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS68333_RigA_20210422T150926';
% expdir0 = '/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS62014_RigD_20210525T133656';
% [~,expname] = fileparts(expdir0);
% expdir = SymbolicCopyExperimentDirectory(expdir0,rootdir);
% results copied to 
% /misc/public/Kristin2Disco/VNC_EXT_VGLUT-GAL4_RigA_20210427T125905/stats.html
% /misc/public/Kristin2Disco/VNC_JRC_SS68333_RigA_20210422T150926/stats.html
% /misc/public/Kristin2Disco/VNC_JRC_SS62014_RigD_20210525T133656/stats.html
for i = 1:numel(explist)
expdir = explist{i};
% FlyDiscoRegisterTrx(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoDectectIndicatorLedOnOff(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoClassifySex(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoComputePerFrameFeatures(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoComputePerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'dorecompute',true,'docomputehists',true,'debugplot',3); % 
% FlyDiscoPlotPerFrameStats2(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'debug',false);
FlyDiscoPlotPerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'debug',false,'makestimvideos',0,'plotstim',1,'plothist',1,'plotflies',true,'plotstimtrajs',1);
end
%%

sd = load(fullfile(expdir,'stats_perframe.mat'));