if ispc,
  rootreaddir = 'E:\Data\FlyBowl\polycarbonate_scratched';
  rootwritedir = 'E:\Data\FlyBowl\polycarbonate_scratched';
  expdir = 'DL-wildtype_TrpA_Rig1Plate01BowlC_20101112T152624';
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  nframes = 9053;
else
  
end
nr = 1024;
nc = 1024;
ncolors = 1;

%%
vidinfo = struct('nr',nr,'nc',nc,'ncolors',ncolors,'nframes',nframes);
obj = Trx('rootreaddir',rootreaddir,'rootwritedir',rootwritedir);
obj.AddExpDir(expdir,vidinfo);