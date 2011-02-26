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
[expdirs_done,expdir_reads_done] = ...
  getExperimentDirs('protocol','20110211','subreadfiles',{'analysis_plots/hist_theta_mm.png'});

expdir_reads = setdiff(expdir_reads,expdir_reads_done);

fid = fopen('expdirs20110215.txt','w');
for i = numel(expdir_reads):-1:1,
  fprintf(fid,'%s\n',expdir_reads{i});
end
fclose(fid);

%% all experiments with stats_perframe data
[expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
  getExperimentDirs('protocol','20110211','subreadfiles',{'stats_perframe.mat'});

fid = fopen('expdirs20110217plot.txt','w');
for i = numel(expdir_reads):-1:1,
  fprintf(fid,'%s\n',expdir_reads{i});
end
fclose(fid);

%% all experiments with hists, stats but no plots

[~,expdirs_ready] = ...
  getExperimentDirs('protocol',analysis_protocol,'subreadfiles',{'hist_perframe.mat'});
[~,expdirs_done] = ...
  getExperimentDirs('protocol',analysis_protocol,'subreadfiles',{'analysis_plots/hist_velmag.png'});
expdirs = setdiff(expdirs_ready,expdirs_done);

fid = fopen('expdirs20110220plot.txt','w');
for i = numel(expdirs):-1:1,
  fprintf(fid,'%s\n',expdirs{i});
end
fclose(fid);

%%

linenames = {'GMR_15D07*','GMR_15H01*','GMR_21C09*','GMR_14G05*','GMR_16C05*',...
  'GMR_16E09*','GMR_16F09*'};
fid = fopen('expdirs20110215avi.txt','w');
for i = 1:numel(linenames),

  [expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
    getExperimentDirs('protocol','20110211','subreadfiles',{'analysis_plots/hist_theta_mm.png'},'linename',linenames{i});
  for j = numel(expdir_reads):-1:1,
    fprintf(fid,'%s\n',expdir_reads{j});
  end
  
end

fclose(fid);