%% set up paths

if ispc,

  addpath E:\Code\SAGE\MATLABInterface\Trunk;
  addpath(genpath('E:\Code\cross_assay\matlab'));
  addpath(genpath('E:\Code\box\PostAnalysis'));
  addpath E:\Code\JCtrax\misc;
  rmSvnPath;
  rootdir = 'E:\Data\FlyBowl\bowl_data';
  
else
  
  addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
  addpath(genpath('/groups/branson/bransonlab/projects/olympiad/cross_assay/matlab'));
  addpath(genpath('/groups/branson/bransonlab/projects/olympiad/box/PostAnalysis'));
  rmSvnPath;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  
end

%% try the data selector

dsf = SAGE.Lab('olympiad').assay('bowl');
dataPull = @(x) pullBowlData(x,'dataset','score','rootdir',rootdir);
selector = OlyDat.DataSelector(dsf,dataPull);

%% load in data chosen with data selector

load('datacache20110420.mat');

%% try the browser
browser = startBowlBrowser(data);