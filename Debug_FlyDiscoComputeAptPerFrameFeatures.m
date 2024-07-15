%% setpath
modpath

%%

settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis';
analysis_protocol = 'settings-internal/20240521_flybubble_LED_VNC3';
params = {'settingsdir',settingsdir,...
    'analysis_protocol',analysis_protocol};

%% expdir
% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240521_testAPTfeatures/VNC2_YNA_K_162984_RigB_20220831T124607';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240605_sanitychecksbouts/VNC3_YNA_K_162984_RigC_20240509T104612';
% 
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_YNA_K_162984_RigD_20210922T133035'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_YNA_K_162984_RigC_20210922T132938'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_YNA_K_162984_RigB_20210922T132823'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_YNA_K_162984_RigA_20210922T132731'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_YNA_K_162984_RigD_20210921T140825'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_YNA_K_162984_RigC_20210921T140740'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_YNA_K_162984_RigB_20210921T140516'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_YNA_K_162984_RigA_20210921T140433'};

%don't walk - need to fix code for them! 
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_JRC_SS68333_RigA_20220906T095615'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_JRC_SS68333_RigA_20220907T121156'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_JRC_SS68333_RigB_20220906T095716'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_JRC_SS68333_RigB_20220907T121238'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_JRC_SS68333_RigC_20220906T095815'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_JRC_SS68333_RigC_20220907T121338'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_JRC_SS68333_RigD_20220906T095905'};
explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS71988_RigC_20210914T133701'
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS71988_RigB_20210914T143211'
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS71988_RigA_20210914T143410'
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS71988_RigD_20210915T155000'
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS71988_RigC_20210915T154857'
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS71988_RigB_20210915T154750'
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS71988_RigA_20210915T154649'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS74168_RigD_20210922T151551'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS74168_RigC_20210922T151455'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS74168_RigB_20210922T151338'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS74168_RigA_20210922T151245'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS74168_RigD_20210921T145837'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS74168_RigC_20210921T145801'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS74168_RigB_20210921T145603'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS74168_RigA_20210921T145515'};

% explist =  {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240521_testAPTfeatures/VNC2_YNA_K_162984_RigB_20220831T124607'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240521_testAPTfeatures/VNC2_YNA_K_162984_RigB_20220831T124607_noaptfeatures'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC_JRC_SS71988_RigC_20210914T133701'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_JRC_SS68333_RigA_20220907T121156'};
%% run code
for i = 1:numel(explist)
    expdir = explist{i}
    FlyDiscoComputeAptPerFrameFeatures(expdir,params{:})
end
