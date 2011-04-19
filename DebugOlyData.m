%% set up paths

addpath E:\Code\SAGE\MATLABInterface\Trunk;
addpath(genpath('E:\Code\cross_assay\matlab'));
addpath(genpath('E:\Code\box\PostAnalysis'));
addpath E:\Code\JCtrax\misc;
rmSvnPath;

%% try the data selector

dsf = SAGE.Lab('olympiad').assay('bowl');
dataPull = @(x) pullBowlData(x,'dataset','score','rootdir','E:\Data\FlyBowl\bowl_data');
selector = OlyDat.DataSelector(dsf,dataPull);
