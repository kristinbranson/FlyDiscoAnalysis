% test whole pipeline

%% set up paths


[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',
    
    addpath E:\Code\JCtrax\misc;
    addpath E:\Code\JCtrax\filehandling;
    addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
    addpath E:\Code\hmm;
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  
  case 'bransonk-lw2',

    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath 'C:\Code\SAGE\MATLABInterface\Trunk';
    addpath C:\Code\hmm;
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
  case 'bransonk-desktop',
    
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;
    addpath /groups/branson/home/bransonk/tracking/code/lds/hmm;
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
end

%% parameters

% expdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_housing_pBDPGAL4U_20111216/results/pBDPGAL4U_UAS_IVS_myr_3_0009_Rig1Plate15BowlA_20120113T142907';
expdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_mating_galit_CS_20120211/results/EXT_CantonS_1220002_None_Rig1Plate15BowlA_20120211T123840';
params = {};

%% run

[success,msgs,stage] = FlyBowlAnalysisPipeline(expdir,params{:});
