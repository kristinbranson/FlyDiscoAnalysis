% make softlinks in the flydisco/to-process folder
%% %%%%% need to run as bransonlab %%%%%%
%% make list
%Katie
% load /groups/branson/home/robiea/Projects_data/Katie/TrpAforReprocessing.mat
% explist = {metadata.file_system_path};
%VNC2 20220609 - wrong gender
% expnamelist = {'VNC2_JRC_SS83407_RigB_20220419T092554',...
% 'VNC2_JRC_SS83407_RigC_20220419T092709',...
% 'VNC2_JRC_SS83407_RigD_20220419T092834',...
% 'VNC2_JRC_SS80084_RigA_20220419T093328',...
% 'VNC2_JRC_SS80084_RigB_20220419T093446',...
% 'VNC2_JRC_SS80084_RigC_20220419T093600',...
% 'VNC2_JRC_SS80084_RigD_20220419T093720',...
% 'VNC2_JRC_SS91980_RigA_20220419T094055',...
% 'VNC2_JRC_SS91980_RigB_20220419T094244',...
% 'VNC2_JRC_SS91980_RigC_20220419T094458',...
% 'VNC2_JRC_SS91980_RigD_20220419T094624',...
% 'VNC2_EXT_VGLUT-GAL4_RigA_20220419T094935',...
% 'VNC2_EXT_VGLUT-GAL4_RigB_20220419T095045',...
% 'VNC2_EXT_VGLUT-GAL4_RigC_20220419T095222',...
% 'VNC2_EXT_VGLUT-GAL4_RigD_20220419T095418',...
% 'VNC2_JRC_SS90380_RigB_20220419T095853',...
% 'VNC2_JRC_SS90380_RigC_20220419T100024',...
% 'VNC2_JRC_SS90380_RigD_20220419T100216',...
% 'VNC2_JRC_SS90380_RigA_20220419T100514',...
% 'VNC2_JRC_SS91891_RigA_20220419T102044',...
% 'VNC2_JRC_SS92103_RigA_20220419T102951'};
% rootdatadir = '/groups/branson/bransonlab/flydisco_data';
% for i = 1:numel(expnamelist)
% explist{i} = fullfile(rootdatadir,expnamelist{i});
% end
%reprocess AmpRec with settingdirectory erros 20220824
% explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/explist_AmpRec_reprocess','%s');
% explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/explist_AmpRec_reprocess_20220825','%s');
% explist = {'/groups/branson/bransonlab/flydisco_data/VNC2_JRC_SS74168_RigD_20230620T112923'
% '/groups/branson/bransonlab/flydisco_data/VNC2_JRC_SS74168_RigC_20230620T112856'
% '/groups/branson/bransonlab/flydisco_data/VNC2_JRC_SS73771_RigD_20230620T122726'
% '/groups/branson/bransonlab/flydisco_data/VNC2_JRC_SS73787_RigD_20230620T111749'};
% explist = {'/groups/branson/bransonlab/flydisco_data/LPC1CsChr_JHS_K_85321_RigB_20231116T132617'
% '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_BJD_SS02407_RigA_20231116T132456'
% '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_GMR_SS02575_RigD_20231116T132939'
% '/groups/branson/bransonlab/flydisco_data/LPC1ShiKir_BJD_SS02407_RigA_20231116T150618'
% '/groups/branson/bransonlab/flydisco_data/LPC1ShiKir_BJD_SS02407_RigC_20231116T151050'
% '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_GMR_SS02575_RigC_20231116T132833'
% '/groups/branson/bransonlab/flydisco_data/LPC1ShiKir_BJD_SS02407_RigD_20231116T151129'
% '/groups/branson/bransonlab/flydisco_data/LPC1ShiKir_BJD_SS02407_RigB_20231116T150803'};
% explist = {'/groups/branson/bransonlab/flydisco_data/LPC1CsChr_JHS_K_85321_RigB_20231116T132617'
% '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_GMR_SS02575_RigD_20231116T132939'
% '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_BJD_SS02407_RigA_20231116T132456'
% '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_GMR_SS02575_RigC_20231116T132833'};

%LPC1 fixed screen_types 
% explist = {'/groups/branson/bransonlab/flydisco_data/LPC1CsChr_GMR_SS02575_RigA_20231120T155527'
% '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_BJD_SS02407_RigD_20231120T153809'
% '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_BJD_SS02407_RigC_20231120T153725'
% '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_JHS_K_85321_RigB_20231120T153618'
% '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_JHS_K_85321_RigA_20231120T153353'};

% % LCP1 fixed missing lbl file
% explist = {'    /groups/branson/bransonlab/flydisco_data/LPC1CsChr_BJD_SS02407_RigA_20231211T140307'
%     '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_JHS_K_85321_RigA_20231211T142439'
%     '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_BJD_SS02407_RigB_20231211T140415'
%     '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_BJD_SS02407_RigC_20231211T140526'
%     '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_BJD_SS02407_RigD_20231211T140651'
%     '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_GMR_SS02575_RigD_20231211T142705'
%     '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_JHS_K_85321_RigB_20231211T142434'
%     '/groups/branson/bransonlab/flydisco_data/LPC1CsChr_JHS_K_85321_RigC_20231211T142559'};

% fixing Kai's DN setting dirs
explist = {'/groups/branson/bransonlab/flydisco_data/DNcontrol_GMR_SS02575_RigB_20231213T135549'
    '/groups/branson/bransonlab/flydisco_data/DNretinal_SS77821_RigA_20231213T141155'};
%%
toprocessdir =  '/groups/branson/bransonlab/flydisco_data/to-process';
%ln -s source_file_or_directory_name  softlink_name
for i = 1:numel(explist)
    expdir = explist{i};
    [~,expname] = fileparts(expdir);
    outdir = fullfile(toprocessdir,expname);
    cmd = sprintf('ln -s %s %s', expdir, outdir);
    system(cmd)
end