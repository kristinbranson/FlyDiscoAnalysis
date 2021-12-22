%Debug_JAABADetect.m 12/16/2016

% %% setup path
% 
% addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/filehandling/;
% addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/misc;
% addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/perframe/;

settingsdir = '/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis/FlyDiscoAnalysis-new/settings';
% change to mine

%% experiment directories
% expdir = '/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis_ExptsToAnalyze_Test/aIPg/20210909T080233_rig1_flyBowl2__20XUASCsChrimsonattp18_SS36564_KS_redonly3times10AND30_080521';
% change to exptdir
expdir = {'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T075256_rig1_flyBowl2__CantonSstvd_CantonS_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T075256_rig1_flyBowl3__OreRstvd_OreR_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T081124_rig1_flyBowl2__OreR_OreR_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T081124_rig1_flyBowl3__NorpA_norpAempty_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T083000_rig1_flyBowl2__CantonS_CantonS_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T083000_rig1_flyBowl3__Piezo_Piezo_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T084743_rig1_flyBowl2__OreRstvd_OreR_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T084743_rig1_flyBowl3__CantonSstvd_CantonS_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T090619_rig1_flyBowl2__norpAempty_NorpA_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T090619_rig1_flyBowl3__OreR_OreR_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T092442_rig1_flyBowl2__Piezo_Piezo_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T092442_rig1_flyBowl3__Orco2_Orco2_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T094405_rig1_flyBowl2__Orco2_Orco2_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/KS_DataFromRigs/FlyBowlRGB/FlyBowl_NonOpto/20211208/20211208T094405_rig1_flyBowl3__CantonS_CantonS_20211001_nonopto'};

%% parameters
analysis_protocol = '20211112_FlyBowlRGBVision_addJAABA';
% change to yours
% analysis_protocol = 'current_non_olympiad_rubin_grooming';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};
%%

datalocparamsfilestr = 'dataloc_params.txt';
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
jaabaclassifierparamsfilestrs = fullfile(settingsdir,analysis_protocol,dataloc_params.jaabaclassifierparamsfilestrs);

jabfiles = textread(jaabaclassifierparamsfilestrs,'%s');

pwdprev = pwd;
jaabadir = fileparts(which('JAABADetect'));
cd(jaabadir);

JAABADetect(expdir,'jabfiles',jabfiles,'forcecompute',false);

cd(pwdprev);