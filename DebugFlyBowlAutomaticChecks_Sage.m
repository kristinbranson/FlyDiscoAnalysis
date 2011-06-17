%% set up path


[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',
    
    addpath E:\Code\JCtrax\misc;
    addpath E:\Code\JCtrax\filehandling;
    addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  
  case 'bransonk-lw2',

    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
  case 'bransonk-desktop',
    
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
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

%% get all experiments

data = SAGEListBowlExperiments('checkflags',false,'rootdir',rootdatadir,'daterange',{'20110201T000000'});

%%

successes = false(1,numel(data));
msgs = cell(1,numel(data));

%%

for i = 1:numel(data),
  fprintf('%d: %s...\n',i,data(i).experiment_name);
  expdir = data(i).file_system_path;
  [successes(i),msgs{i}] = FlyBowlAutomaticChecks_Sage(expdir,'dosave',false);
end

%% print to file

fid = fopen('AutomaticChecksSageResults20110616.txt','w');
for i = 1:numel(data),
  if successes(i),
    fprintf(fid,'%s,P\n',data(i).experiment_name);
  else
    fprintf(fid,'%s,F\n',data(i).experiment_name);
  end
end
fclose(fid);
