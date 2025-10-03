modpath
% 
% % %run code with KB's March 2024 updates
% % % explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240315_statplotting_KBupdates/VNC2_EXT_VGLUT-GAL4_RigB_20231114T124058'};
% % % explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240315_statplotting_KBupdates/VNC2_JRC_SS39036_RigB_20231024T114951'};
% % explist= {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240315_statplotting_KBupdates/LPC1CsChr_GMR_SS02575_RigA_20231120T155527'};
% settingsdir = 'settings-internal';
% 
% % % analysis_protocol = '20230907_flybubble_LED_VNC2_hack';
% 
% % LPC1
% % explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240315_statplotting_KBupdates/LPC1CsChr_BJD_SS02407_RigC_20231211T140526'};
% % analysis_protocol = '20240402_flybubble_LED_LPC1_CsChr';
% 
% % % VNC
% % % % explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240315_statplotting_KBupdates/VNC2_EXT_VGLUT-GAL4_RigB_20230516T125145'};
% % explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240315_statplotting_KBupdates/VNC2_EXT_VGLUT-GAL4_RigB_20230705T111622'};
% % explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20250328_statsperframe/VNC2_EXT_VGLUT-GAL4_RigA_20230725T123307'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20250328_statsperframe/VNC2_EXT_VGLUT-GAL4_RigB_20230705T111622'};
% analysis_protocol = '20240402_flybubble_LED_VNC2';

% multi indicators
settingsdir = '/groups/branson/home/robiea/Projects_data/Adriane/multipleindicators/alice_testing';
analysis_protocol = '20250306_MaleFemaleGtACRRig3';
explist = {'/groups/branson/home/robiea/Projects_data/Adriane/multipleindicators/alice_testing/MaleFemaleGtACRRig3_multibubble__P1GtACR_SS102716_UASGtACR_20250314T111250'};

% % explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240404_statsplotting_Katie/20230719T074540_rig1_flyBubble1_SS36564_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR'};
% % analysis_protocol ='20240329_BubblefRGB_GtAgStats_AR';

% for i = 1:numel(explist)
% expdir = explist{i};
% % FlyDiscoRegisterTrx(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% % FlyDiscoDetectIndicatorLedOnOff(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% % FlyDiscoClassifySex(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% % FlyDiscoComputePerFrameFeatures(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoComputePerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'dorecompute',false,'docomputehists',true,'debugplot',false); % 
% FlyDiscoPlotPerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'debug',false,'makestimvideos',0,'plotstim',true,'plothist',false,'plotflies',false,'plotstimtrajs',false);
% end

% FlyDiscoComputePerFrameStats
% dorecompute <true or false>  force recompute or only recompute stats if they don't exist
%debugplot,false or n plots - plots the firs n stat_perframefeatures to
%visualizeanalyzed frames and feature values 

% FlyDiscoPlotPerFrameStats
% debug <true or false> true displays plots but does not save them, false
% saves plots but not display them
% TO DO
% makestimvideos
% plotstim
% plothist
% plotflies
% plotstimtrajs

%run everything
for i = 1:numel(explist)
expdir = explist{i};
% FlyDiscoRegisterTrx(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoDetectIndicatorLedOnOff(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoClassifySex(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoComputePerFrameFeatures(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoComputePerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'dorecompute',true,'docomputehists',true,'debugplot',2);  
FlyDiscoPlotPerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'debug',false,'plotstats',1,'makestimvideos',1,'plotstim',1,'plothist',1,'plotflies',1,'plotstimtrajs',1,'plottimeseries',1);
end
%%

% sd = load(fullfile(expdir,'stats_perframe.mat'));