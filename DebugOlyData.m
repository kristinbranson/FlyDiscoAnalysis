%% set up paths

if ispc,

  addpath E:\Code\SAGE\MATLABInterface\Trunk;
  addpath(genpath('E:\Code\cross_assay\matlab'));
  addpath(genpath('E:\Code\box\PostAnalysis'));
  addpath E:\Code\JCtrax\misc;
  rmSvnPath;
  
else
  
  addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
  addpath(genpath('/groups/branson/bransonlab/projects/olympiad/cross_assay/matlab'));
  addpath(genpath('/groups/branson/bransonlab/projects/olympiad/box/PostAnalysis'));
  rmSvnPath;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;

end

%% try the data selector

dsf = SAGE.Lab('olympiad').assay('bowl');
dataPull = @(x) pullBowlData(x,'dataset','score','rootdir','E:\Data\FlyBowl\bowl_data');
selector = OlyDat.DataSelector(dsf,dataPull);
