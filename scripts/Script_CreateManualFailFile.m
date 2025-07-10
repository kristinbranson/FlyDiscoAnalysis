%% write manual fail files 
%%%must run as bransonlab for VNC screen data%%%
%% set path
% modpath

%% exp list and failure reason (keep failure reason list)
% csv file for manual fails 
failtablefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/manualfailures_log.csv';
explist = [];
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20221107_testingmetadatafixes/VNC2_YNA_K_162984_RigD_20221004T123313'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20220913_testingJAABA/VNC_JRC_SS46677_RigC_20210602T154857'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20220913_testingJAABA/VNC2_JRC_SS83065_RigA_20220531T104126'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20220913_testingJAABA/VNC2_JRC_SS90562_RigD_20220906T124140'};
% failure_category  = 'testing code';
% 

% rootdatadir 
rootdatadir = '/groups/branson/bransonlab/flydisco_data';
% from reviewing videos with >=12 trajectories (2023 batch)
% % crub in bubble 
% expdirlist = {'VNC2_JRC_SS78430_RigA_20230627T123728'
% 'VNC2_JRC_SS74350_RigB_20230627T122655'
% 'VNC2_JRC_SS50325_RigD_20230830T115927'}
% failure_category = 'crud in bubble';

% % tracked fly outside the bubble
% expdirlist = {'VNC2_EXT_VGLUT-GAL4_RigB_20230621T114253'};
% failure_category =  'tracked fly outside bubble';

% % tracked dead or damaged fly
% expdirlist = {'VNC2_EXT_VGLUT-GAL4_RigC_20230725T123503'
% 'VNC2_EXT_VGLUT-GAL4_RigA_20230627T105833'
% 'VNC2_JRC_SS77231_RigA_20230627T114356'};
% failure_category = 'tracked dead or damaged fly';

% prescreen practicing for rig 4/19/2023 to 5/03/2023
% expdirlist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirlist_prescreening_practice. txt','%s');
% failure_category = 'testing';

% WL 43
% expdirlist = {'VNC2_EXT_VGLUT-GAL4_RigA_20231205T115555'};
% failure_category =  'crud in bubble';
% % WL # 38
% expdirlist = {'VNC2_JRC_SS76743_RigB_20230926T111421'};
% failure_category =  'crud in bubble';
% % WL # 34
% expdirlist = {'VNC2_EXT_VGLUT-GAL4_RigC_20230830T121051'};
% failure_category =  'bumped rig';
% % % WL 33
% expdirlist = {'VNC2_YNA_K_162984_RigA_20230725T121228'};
% failure_category =  'crud in bubble';
% expdirlist = {'VNC2_JRC_SS92097_RigA_20230725T140635'
% 'VNC2_JRC_SS92097_RigB_20230725T140655'
% 'VNC2_JRC_SS92097_RigC_20230725T140808'
% 'VNC2_JRC_SS92097_RigD_20230725T140824'}
% failure_category = 'late runs';

% WL #32
% expdirlist = {'VNC2_JRC_SS90441_RigA_20230718T140114'
% 'VNC2_JRC_SS90441_RigB_20230718T140145'
% 'VNC2_JRC_SS90441_RigC_20230718T140254'
% 'VNC2_JRC_SS90441_RigD_20230718T140333'}
% failure_category = 'late runs';
% % WL 28
% expdirlist = {'VNC2_JRC_SS73860_RigA_20230620T094433'};
% failure_category =  'crud in bubble';

% % WL 26
% expdirlist = {'VNC2_JRC_SS65770_RigC_20230607T095522'};
% failure_category =  'crud in bubble';

% % WL 24
% expdirlist = {'VNC2_JRC_SS83080_RigA_20230516T123901'};
% failure_category =  'crud in bubble';
% % WL 23
% expdirlist = {'VNC2_EXT_VGLUT-GAL4_RigD_20230509T111055'};
% failure_category = 'tracked dead or damaged fly';
% expdirlist = {'VNC2_EXT_VGLUT-GAL4_RigD_20230509T120105'};
% failure_category = 'tracked dead or damaged fly';

% % wk 39 
% expdirlist = {'VNC2_JRC_SS96599_RigA_20220928T125859'
% 'VNC2_JRC_SS96599_RigB_20220928T125932'
% 'VNC2_JRC_SS96599_RigC_20220928T130030'
% 'VNC2_JRC_SS96599_RigD_20220928T130148'};
% 
% failure_category  = 'wrong eye color in cross';
% 
% %wk 36
% expdirlist = {'VNC2_JRC_SS86630_RigC_20220907T122625'};
% failure_category = 'crud in bubble';
% % wk 34
% expdirlist = {'VNC2_JRC_SS92169_RigB_20220803T122924'};
% failure_category = 'crud in bubble';

% %wk 32
% expdirlist = {'VNC2_JRC_SS91177_RigB_20220721T113011'};
% failure_category = 'crud in bubble';

% % effector testing pre week 32
% expdirlist ={'VNC2_EXT_VGLUT-GAL4_RigA_20220713T093804'
% 'VNC2_EXT_VGLUT-GAL4_RigB_20220713T093846'
% 'VNC2_EXT_VGLUT-GAL4_RigC_20220713T094045'};
% failure_category = 'testing';
% 4th effector testing movie 
% expdirlist = {'TrpAFemale4_EXT_VGLUT-GAL4_RigD_20220713T094246'};
% failure_category = 'testing';

% % wk 31 bad effector
% expdirlist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/wk31_manualfail_explist.txt','%s');
% failure_category  = 'bad eye color effector';

% % wk 30 bad effector
% expdirlist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/wk30_manualfail_explist.txt','%s');
% failure_category  = 'bad eye color effector';

% % wk 29 line name errors    
% expdirlist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/wk29_manualfail_explist.txt','%s');
% failure_category  = 'line name error';

% % wk 25 long load time
% expdirlist = {'VNC2_YNA_K_162984_RigB_20220503T111634'};
% failure_category  = 'long load time';

% % wk 24 long load time
% expdirlist = {'VNC2_JRC_SS79392_RigC_20220426T113912'};
% failure_category  = 'long load time';

% % wk 23 long load time
% expdirlist = {'VNC2_JRC_SS90380_RigA_20220419T100514'
%     'VNC2_JRC_SS91891_RigB_20220419T102859'};
% failure_category  = 'long load time';
% expdirlist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/wk23_manualfail_explist.txt','%s');
% failure_category  = 'cold plate too cold';

% % wk 21 long load time
% expdirlist = {'VNC_JRC_SS70427_RigA_20211027T154243'};
% failure_category  = 'long load time';
 
% % wk 17 long load time
% expdirlist = {'VNC_JRC_SS60607_RigB_20210929T150636'};
% failure_category  = 'long load time';

% % wk 16 long load time
% expdirlist = {'VNC_JRC_SS68298_RigB_20210921T143242'};
% failure_category  = 'long load time';

% % wk 12 long load time 
% expdirlist = {'VNC_EXT_VGLUT-GAL4_RigC_20210614T150203'};
% failure_category  = 'long load time';

% % wk 10 crud in bubble
% expdirlist = {'VNC_JRC_SS66958_RigC_20210525T142206'};
% failure_category = 'crud in bubble';

% %wk 9 long load time
% expdirlist = {'VNC_JRC_SS60250_RigA_20210520T133054'};
% failure_category  = 'long load time';

% % wk 7 crud in bubble
% expdirlist = {'VNC_JRC_SS53890_RigB_20210503T141702'};
% failure_category = 'crud in bubble';
% % long load time
% expdirlist = {'VNC_JRC_SS59163_RigB_20210503T134731'};
% failure_category = 'long load time';

% % wk 5 long load time
% expdirlist = {'VNC_JRC_SS51902_RigC_20210420T152938'};
% failure_category = 'long load time';

% % wk 4 long load time
% expdirlist = {'VNC_JRC_SS45584_RigA_20210415T143020'
%     'VNC_JRC_SS44268_RigA_20210414T150229'};
% failure_category = 'long load time';

% % wk 3 long load time
% expdirlist = {'VNC_JRC_SS39036_RigC_20210407T150806'
% 'VNC_JRC_SS38644_RigC_20210405T134827'
% 'VNC_JRC_SS38665_RigA_20210405T155121'
% 'VNC_JRC_SS38649_RigD_20210407T140431'
% };
% failure_category = 'long load time';

% bubble loaded into rig incorrectly
% expdirlist = {'VNC_JRC_SS38679_RigA_20210405T160207'};
% failure_category = 'bubble placement error';

% % wk 2 long load time
% expdirlist = {'VNC_JRC_SS38592_RigD_20210330T142057'
% 'VNC_JHS_K_85321_RigD_20210331T150553'};
% failure_category = 'long load time';
% % wk 2 more
% expdirlist = {'VNC_JRC_SS37641_RigC_20210401T140146'
%     'VNC_JRC_SS38337_RigD_20210401T151839' % already deleted for aborted? 
%     'VNC_EXT_VGLUT-GAL4_RigD_20210330T135042'};
% failure_category = 'long load time';
% %wk 2 more
% expdirlist = {'VNC_JRC_SS36165_RigD_20210330T133117'}
% failure_category = 'tracked dead or damaged fly';
% % wk 1 fail starvation too long 
% expdirlist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/wk1_manualfail_explist.txt','%s');
% failure_category = 'starvation too long';

% % found by Kristin during APT tracking
% expdirlist = {'VNC_JRC_SS33416_RigD_20210601T141610'};
% failure_category = 'bad video';

% expdirlist = {'VNC2_JRC_SS94945_RigB_20220913T122904'};
% failure_category = 'crud in bubble';

% from review of expdirs with trajnum >=12
% % bubble occluded
% expdirlist = {'VNC2_JRC_SS77680_RigA_20220524T111405'};
% failure_category = 'bubble occluded';

% % crud in bubble
% expdirlist = {'VNC2_JRC_SS85803_RigD_20220511T105339'};
% failure_category = 'crud in bubble';

% % dead or damaged flies
% expdirlist = {'VNC_JHS_K_85321_RigC_20210405T140545'
% 'VNC_YNA_K_162984_RigD_20210406T150930'
% 'VNC_EXT_VGLUT-GAL4_RigA_20210413T145019'
% 'VNC_JRC_SS44261_RigB_20210412T151441'
% 'VNC_JRC_SS37648_RigD_20210401T141805'
% 'VNC_EXT_VGLUT-GAL4_RigB_20210628T152220'
% 'VNC_JRC_SS40619_RigD_20210406T134659'
% 'VNC2_JRC_SS74638_RigD_20220726T112118'
% 'VNC2_YNA_K_162984_RigA_20221005T120352'
% 'VNC2_JRC_SS68904_RigC_20220913T103702'
% 'VNC2_EXT_VGLUT-GAL4_RigA_20220830T101614'};
% failure_category = 'tracked dead or damaged fly';

% % fly on outside of bubble
% expdirlist = {'VNC_JRC_SS34580_RigB_20210504T141103'
% 'VNC_JRC_SS50326_RigA_20210427T131657'
% 'VNC_EXT_VGLUT-GAL4_RigC_20210601T143220'
% 'VNC_JRC_SS38601_RigD_20210601T144117'};
% failure_category = 'tracked fly outside bubble';

% fail WL 51 - temperature issue July 2024
% expdirlist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/julytempissue_20240709and10.txt','%s');
% failure_category = 'temperature ctrl failed';

% WL 45
% expdirlist = {'VNC3_JRC_SS34708_RigD_20240515T122741'};
% failure_category = 'crud in bubble';
%WL 46;
% expdirlist = {'VNC3_JRC_SS46737_RigA_20240521T113401'};
% failure_category = 'crud in bubble';
% WL 47
% expdirlist = {'VNC3_JRC_SS57959_RigB_20240528T111825'
%     'VNC3_JRC_SS52765_RigB_20240528T121058'
%     'VNC3_JRC_SS59153_RigB_20240528T125232'};
% 
% failure_category = 'crud in bubble';
% WL 50
% expdirlist = {'VNC3_JRC_SS37671_RigB_20240702T131521'};
% failure_category  = 'long load time';
% expdirlist = {'VNC3_JRC_SS38013_RigA_20240702T111149'
%     'VNC3_JRC_SS38013_RigC_20240702T11133'
%     'VNC3_JRC_SS38631_RigA_20240702T112317'};
% expdirlist = {'VNC3_JRC_SS46233_RigD_20240702T120903'};
% expdirlist = {'VNC3_JRC_SS38013_RigC_20240702T111337'}
% failure_category = 'crud in bubble';

%WL 52
% expdirlist = {'VNC3_JRC_SS59230_RigC_20240723T104809'
% 'VNC3_JRC_SS65765_RigA_20240723T110819'
% 'VNC3_JRC_SS65765_RigD_20240723T111041'
% 'VNC3_JRC_SS53050_RigA_20240723T112105'
% 'VNC3_JRC_SS59225_RigA_20240723T114406'
% 'VNC3_JRC_SS61036_RigA_20240723T115702' };
% failure_category = 'crud in bubble';
% %WL 53
% expdirlist = {'VNC3_EXT_VGLUT-GAL4_RigC_20240730T110155'
% 'VNC3_EXT_VGLUT-GAL4_RigD_20240730T110310'};
% failure_category = 'tracked dead or damaged fly';
% WL 54
% expdirlist = {'VNC3_YNA_K_162984_RigD_20240806T115455'};
% failure_category = 'crud in bubble';
% WL 55
% expdirlist = {'VNC3_JRC_SS105097_RigA_20240813T115528'};
% failure_category = 'crud in bubble';
% WL 56
% expdirlist = {'VNC3_JRC_SS102245_RigA_20240827T122718'
%     'VNC3_JRC_SS96563_RigC_20240827T111512'};
% failure_category = 'crud in bubble';
% WL 59
% expdirlist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/septtempissue_20240910.txt','%s');
% failure_category = 'temperature ctrl failed';
% WL 61
% expdirlist = {'VNC3_JRC_SS50974_RigA_20240924T103307'};
% failure_category = 'crud in bubble';
% WL 62
% expdirlist = {'VNC3_JRC_SS100065_RigB_20241008T112542'
% 'VNC3_EXT_VGLUT-GAL4_RigA_20241009T105712'
% 'VNC3_JRC_SS100065_RigB_20241009T121453'};
% failure_category = 'crud in bubble';

% WL 63
% expdirlist = {'VNC3_JRC_SS102259_RigA_20241016T104603'
% 'VNC3_JRC_SS102311_RigD_20241016T115405'
% 'VNC3_JRC_SS48644_RigC_20241017T114320'};
% failure_category = 'crud in bubble';

% reviewing high trajectory number
% expdirlist = {'VNC3_JRC_SS46706_RigD_20240924T121956'};
% failure_category = 'tracked dead or damaged fly';
% expdirlist = {'VNC3_JRC_SS97456_RigA_20240807T114733'};
% failure_category = 'crud in bubble';
expdirlist = {'VNC3_JRC_SS100086_RigA_20240820T124724'};
failure_category = 'tracked dead or damaged fly';
%%
for i = 1:numel(expdirlist)
explist{i} = fullfile(rootdatadir,expdirlist{i});
end

%% run function

[success] = CreateManualFailFile(explist,failure_category,'Replace',true); 

%% write failures to a csv files as well 

if ~exist(failtablefile,'file')
    fid = fopen(failtablefile,'a');
    fprintf(fid,'expdir\t failure type\t file written\t ');
else
    fid = fopen(failtablefile,'a');
end
for i = 1:numel(explist)
fprintf(fid,'\n%s\t %s\t %d\t ',expdirlist{i},failure_category,success(i));
end

fclose(fid);

