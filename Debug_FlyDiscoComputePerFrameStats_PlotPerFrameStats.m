modpath

%run code with KB's March 2024 updates
explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240315_statplotting_KBupdates/VNC2_EXT_VGLUT-GAL4_RigB_20231114T124058'};
settingsdir = 'settings-internal';
analysis_protocol = '20230907_flybubble_LED_VNC2';



for i = 1:numel(explist)
expdir = explist{i};
% FlyDiscoRegisterTrx(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoDetectIndicatorLedOnOff(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoClassifySex(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% FlyDiscoComputePerFrameFeatures(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
FlyDiscoComputePerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'dorecompute',false,'docomputehists',true,'debugplot',0); % 
FlyDiscoPlotPerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,'debug',false,'makestimvideos',1,'plotstim',true,'plothist',true,'plotflies',true,'plotstimtrajs',true);
end
%%

% sd = load(fullfile(expdir,'stats_perframe.mat'));