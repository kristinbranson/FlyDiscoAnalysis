%% set up paths

addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath(genpath('/groups/branson/bransonlab/projects/olympiad/cross_assay/matlab'));
addpath(genpath('/groups/branson/bransonlab/projects/olympiad/box/PostAnalysis'));
rmSvnPath;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath olydat_browser;

rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
settingsdir = 'settings';

%% set parameters

params = struct;
params.dataset = 'score';
params.data_type = {'sexclassifier_diagnostics_classifier_mu_area_female',...
  'sexclassifier_diagnostics_classifier_mu_area_male'};
params.screen_type = 'primary';
params.daterange = {'20110201T000000'};
params.checkflags = true;
params.removemissingdata = true;
params.rootdir = rootdatadir;

paramscell = struct2paramscell(params);
analysis_protocol = '20111005';

%% pull data

data = SAGEGetBowlData(paramscell{:});
data = PrepareOlyDatData(data);

%%

browser = startBowlBrowser(data,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);