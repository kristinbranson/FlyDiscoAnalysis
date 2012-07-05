%% set up path

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',
    
    addpath E:\Code\JCtrax\misc;
    addpath E:\Code\JCtrax\filehandling;
    addpath E:\Code\FlyBowlDataCapture;
    addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  
  case 'bransonk-lw2',

    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath C:\Code\FlyBowlDataCapture;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
  case 'bransonk-desktop',
    
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlDataCapture;
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath C:\Code\FlyBowlDataCapture;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
end

analysis_protocol = '20120210';
outfilename = 'SAGEtoExpDirConsistencyCheckResults20120430.tsv';

%%

metadata = SAGEListBowlExperiments('removemissingdata',false,'checkflags',false);  

fid = fopen(outfilename,'w');
[~,order] = sort(str2double({metadata.wish_list}));
order = order(end:-1:1);
metadata = metadata(order);

wishlistcurr = '????';
for i = 1:numel(metadata),
  expdir = metadata(i).file_system_path;
  if ~exist(expdir,'dir'),
    if ~strcmp(wishlistcurr,metadata(i).wish_list),
      wishlistcurr = metadata(i).wish_list;
      fprintf(fid,'\n****WISHLIST %s****\n',wishlistcurr);
      fprintf('\n****WISHLIST %s****\n',wishlistcurr);
    end
    fprintf(fid,'%s\n',expdir);
    fprintf('MISSING %s\n',expdir);
  end
end

fclose(fid);

%% set parameters

params = {'settingsdir',settingsdir,'rootdatadir',rootdatadir,...
  'analysis_protocol',analysis_protocol,'outfilename',outfilename,...
  'daterange',{'20110201T000000'}};

%% 

[isconsistent,filesmissing] = SAGEExpDirConsistencyCheckMany(params{:});%,'restartexperiment','GMR_38H09_AE_01_TrpA_Rig1Plate15BowlD_20120105T160100');


% find experiments that are in SAGE that are not on the file system

