%% set up paths
[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',

    addpath E:\Code\SAGE\MATLABInterface\Trunk;
    addpath(genpath('E:\Code\cross_assay\matlab'));
    addpath(genpath('E:\Code\box\PostAnalysis'));
    addpath E:\Code\JCtrax\misc;
    addpath olydat_browser;
    rmSvnPath;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    
  case 'bransonk-lw2',

    addpath C:\Code\SAGE\MATLABInterface\Trunk;
    addpath(genpath('C:\Code\cross_assay\trunk\matlab'));
    addpath(genpath('C:\Code\box\PostAnalysis\trunk\matlab'));
    addpath C:\Code\JCtrax\misc;
    addpath olydat_browser;
    rmSvnPath;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';

  case 'bransonk-desktop',
    
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/cross_assay/matlab'));
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/box/PostAnalysis'));
    rmSvnPath;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath olydat_browser;
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    
  case 'robiea-ww1'
    
    addpath C:\Users\robiea\Documents\Code_versioned\SAGE\MATLABInterface\Trunk;
    addpath(genpath('C:\Users\robiea\Documents\Code_versioned\cross_assay'));
    addpath(genpath('C:\Users\robiea\Documents\Code_versioned\PostAnalysis'));
    addpath C:\Users\robiea\Documents\Code_versioned\JCtrax\misc;
    addpath C:\Users\robiea\Documents\Code_versioned\FlyBowlAnalysis\olydat_browser;
    rmSvnPath;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Users\robiea\Documents\Code_versioned\FlyBowlAnalysis\settings';
    
  case 'robiea-ws'
    addpath '/groups/branson/home/robiea/Code_versioned/SAGE/MATLABInterface/Trunk';
    addpath(genpath('/groups/branson/home/robiea/Code_versioned/OlyDat/matlab'));
    addpath(genpath('/groups/branson/home/robiea/Code_versioned/OlyDat/postanalysis'));
    rmSvnPath
    addpath '/groups/branson/home/robiea/Code_versioned/JCtrax/misc';
    addpath '/groups/branson/home/robiea/Code_versioned/JCtrax/filehandling';
    addpath '/groups/branson/home/robiea/Code_versioned/FlyBowlAnalysis/olydat_browser';
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBowlAnalysis/settings';
    
  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/cross_assay/matlab'));
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/box/PostAnalysis'));
    rmSvnPath;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath olydat_browser;
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    
end

analysis_protocol = '20110804';

%% pull lots of control data

params = struct;
% data set: data -- we want behavior data
params.dataset = 'data';
% only pull data marked as primary screen data
params.screen_type = 'primary';
% only pull control data
params.line_name = 'pBDPGAL4U';
% don't use the automatic flag checking implemented in SAGEGetBowlData,
% we'll do that here if we want to
params.checkflags = false;
% don't remove experiments with missing data, we'll do that manually if we
% decide to. 
params.removemissingdata = false;
% don't pull aborted experiments
params.flag_aborted = 0;
% don't pull experiments with flag_redo set to 1
params.flag_redo = 0;
% restrict to experiments with automated_pf = P
% params.automated_pf = 'P';
% restrict to experiments with manual_pf ~= F
% params.manual_pf = {'P','U'};

paramscell = struct2paramscell(params);

data = {};
day0 = now;
dateformat = 'yyyymmddTHHMMSS';
period = 2;
for daysago = 0:period:200,
  tryn = 1;
  while true,
    success = false;
    try
      fprintf('daysago = %d, try = %d\n',daysago,tryn);
      data{end+1} = SAGEGetBowlData(paramscell{:},'daterange',{datestr(day0-daysago-period,dateformat),datestr(day0-daysago,dateformat)});
      success = true;
    catch ME,
      delete(waitbar(0));
      getReport(ME)
      tryn = tryn + 1;
    end
    if success,
      break;
    end
  end
  save ControlData20110902.mat;
end
