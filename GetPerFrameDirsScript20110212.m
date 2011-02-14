% set up path
if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
end

[expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
  getExperimentDirs('protocol','20110211','subreadfiles',{'perframe/velmag.mat'});
fid = fopen('expdirs20110212.txt','w');
for i = numel(expdir_reads):-1:1,
  fprintf(fid,'%s\n',expdir_reads{i});
end
fclose(fid);