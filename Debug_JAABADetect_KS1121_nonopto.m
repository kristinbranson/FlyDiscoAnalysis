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
expdir = {'/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis_ExptsToAnalyze_Test/2021-11-05/20211105T080526_rig1_flyBowl2__NorpA_norpAempty_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis_ExptsToAnalyze_Test/2021-11-05/20211105T082605_rig1_flyBowl1__CantonS_CantonS_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis_ExptsToAnalyze_Test/2021-11-05/20211105T082605_rig1_flyBowl2__Piezo_Piezo_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis_ExptsToAnalyze_Test/2021-11-05/20211105T084154_rig1_flyBowl1__CantonS_CantonS_20211001_nonopto',...
'/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis_ExptsToAnalyze_Test/2021-11-05/20211105T085601_rig1_flyBowl1__Piezo_Piezo_20211001_nonopto'};

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